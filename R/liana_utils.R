

#' Function to translate Human genes to organism X orthologues
#'
#' @param op_resource resource in the format of OmniPath/LIANA
#'
#' @param symbols_dict dictionary (named list) with human genesymbols and ortholog
#'  in a second species
#'
#' @param .default_fun default function to be applied to any missing/unmapped
#' genesymbols (`str_to_title` by default)
#'
#' @export
generate_orthologs <- function(op_resource,
                               symbols_dict,
                               .default_fun = "str_to_title"){
    # decomplexify
    op_resource_decomplex <- op_resource %>%
        decomplexify()

    # translate subunits
    translated_subunits <- op_resource_decomplex %>%
        mutate(across(ends_with("genesymbol"),
                      ~recode.character2(.x, !!!symbols_dict,
                                         .default_fun = .default_fun)))

    # Generate Dictionaries for complexes
    target_complex_dict <- .generate_complex_dict(translated_subunits,
                                                  entity="target") %>%
        deframe()
    source_complex_dict <- .generate_complex_dict(translated_subunits,
                                                  entity="source") %>%
        deframe()

    # Bind all dictionaries
    dict <- pmap( # append multiple lists
        list(
            list(
                symbols_dict,
                target_complex_dict,
                source_complex_dict
            )
        ), c) %>%
        flatten %>%
        flatten()

    # get orthologous resource
    op_ortholog <- op_resource %>%
        mutate(across(ends_with("genesymbol"),
                      ~recode.character2(.x, !!!dict,
                                         .default_fun = .default_fun)))

    return(op_ortholog)
}


#' Modified `dplyr::recode` function
#'
#' @inheritParams dplyr::recode
#'
#' @param .default_fun Function to modify the default (mismatched) strings
#' (`str_to_title` by default - i.e. for mice)
#'
#' @details enables to modify unmatched genesymbols
#' @import stringr
#'
#' @noRd
recode.character2 <- function(.x,
                              ...,
                              .default = NULL,
                              .missing = NULL,
                              .default_fun) {
    .x <- as.character(.x)
    values <- rlang::list2(...)
    if (!all(rlang::have_name(values))) {
        bad <- which(!rlang::have_name(values)) + 1
        msg <- glue::glue("{dplyr:::fmt_pos_args(bad)} must be named.")
        rlang::abort(msg)
    }

    n <- length(.x)
    template <- dplyr:::find_template(values, .default, .missing)
    out <- template[rep(NA_integer_, n)]
    replaced <- rep(FALSE, n)

    for (nm in names(values)) {
        out <- dplyr:::replace_with(out, .x == nm, values[[nm]], paste0("`", nm, "`"))
        replaced[.x == nm] <- TRUE
    }

    .default <- dplyr:::validate_recode_default(.default, .x, out, replaced)
    out <- dplyr:::replace_with(out, !replaced & !is.na(.x), exec(.default_fun, .default), "`.default`")
    out <- dplyr:::replace_with(out, is.na(.x), .missing, "`.missing`")
    out
}




#' Helper function to generate a dictionary also for the complexes
#'
#' @param translated_subunits decomplexified op_resource with already recoded
#' subunits
#'
#' @param entity column of interest (target or source)
#'
#' @noRd
.generate_complex_dict <- function(translated_subunits, entity = "target"){
    genesymbol_complex <- str_glue("{entity}_genesymbol_complex")
    genesymbol <- str_glue("{entity}_genesymbol")

    filter(translated_subunits, str_detect(.data[[genesymbol_complex]], "_")) %>%
        dplyr::select(.data[[genesymbol_complex]], .data[[genesymbol]]) %>%
        distinct() %>%
        group_by(.data[[genesymbol_complex]]) %>%
        group_nest(keep = FALSE, .key = 'subunits') %>%
        mutate(translated_complex = map(subunits, function(sub){
            glue::glue_collapse(sub[[genesymbol]], sep = "_") %>%
                as.character()
        })) %>%
        dplyr::select(.data[[genesymbol_complex]], translated_complex) %>%
        unnest(translated_complex)

}

