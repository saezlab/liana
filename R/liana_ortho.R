
#' Function to generate a homologous OmniPath resource
#'
#' @param op_resource a resource in the format of OmniPath/LIANA
#'
#' @param target_organism `ncbi_taxid` or `name` of the target organism.
#' See `show_homologene` for available organisms via OmnipathR's `HomoloGene`
#'
#' @param max_homologs Determines the max number of homologs to be translated.
#' Certain genes will have multiple homolog matches, with some having also
#' certain isoforms considered. To exclude cases in which the number of
#' matched homologs is too high, one can adjust the homologs parameter.
#' Setting this to `1` would mean that one-to-many homolog matches are discarded
#'
#' @param .missing_fun approach to handle missing interactions. By default
#' set to `NULL` which would mean that any interactions without a homology
#' match will be filtered. This can be set to e.g. `str_to_title` when working
#' with murine symbols. Then if a gene has no matched homolog, instead of
#' discarding it, the `.missing_fun` will be used to format the name from human.
#' Hence, increasing the number of matches, but likely introducing some
#' mismatches.
#'
#' @param symbols_dict `NULL` by default, then `homology_dict` is called
#' to generate a dictionary from OmniPathR's homologene resource. Alternatively,
#' one can pass their own symbols_dictionary.
#'
#' @param source name of the source (ligand) column
#'
#' @param target name of the target (receptor) column
#'
#' @param verbose logical for verbosity
#'
#' @param mappings Character vector: control ambiguous mappings.
#'
#' @return a converted ligand-receptor resource
#'
#' @export
generate_homologs <- function(op_resource,
                              target_organism,
                              max_homologs = 5,
                              .missing_fun = NULL,
                              symbols_dict = NULL,
                              columns = c("source_genesymbol",
                                          "target_genesymbol"),
                              verbose = TRUE,
                              mappings = c("1:1", "1:m", "n:1", "n:m")){

    op_resource %<>% mutate(across(all_of(columns),
                                   ~str_replace(., "COMPLEX:", "")))

    # Minimum column set resource
    minres <- op_resource %>% select(!!columns)

    # Get decomplexified resource
    decomp <- decomplexify(minres, columns = columns)

    # Get union of symbols
    entities <- purrr::reduce(map(columns, function(col) decomp[[col]]), union)

    # generate homology geneset
    symbols_dict <-
        homology_dict(
            entities = entities,
            target_organism = target_organism,
            mappings = mappings
        )


    # Remove any missing antities
    if(is.null(.missing_fun)){

        # All missing entities
        missing_entities <- setdiff(entities,
                                    names(symbols_dict))

        liana_message(
            str_glue("Entries without homologs ({length(missing_entities)}):
                     {paste(missing_entities, collapse = '; ')}"),
            verbose = verbose
        )

        # Keep only interactions for which all proteins (incl. subunits) have a matching homologue
        missing <- decomp %>%
            # check if neither is in missing entities
            mutate(lr_present = !if_any(columns, function(x) x %in% missing_entities)) %>%
            group_by(across(all_of(ends_with("complex")))) %>%
            summarise(all_present = mean(lr_present), .groups = "keep") %>%
            # only keep those that are present
            filter(all_present < 1) %>%
            # remove _complex for join
            rename_with(~gsub("_complex", "", .x), ends_with("complex"))

        # Remove interactions without matches
        minres <- anti_join(minres,
                            missing,
                            by = columns)
    }


    # Obtain genes with multiple matches
    entity_2many <- symbols_dict %>%
        enframe(name = "genesymbol_source",
                value = "genesymbol_target") %>%
        group_by(genesymbol_source) %>%
        count(name = "n_match") %>%
        filter(n_match > 1 & n_match <= max_homologs) %>%
        pull(genesymbol_source)

    liana_message(
        stringr::str_glue(
            paste0(
                "One-to-many homolog matches ({length(entity_2many)}): ",
                "{paste(entity_2many, collapse = '; ')}"
            )
        ),
        verbose = verbose
    )

    ### IF NOT many2many .handle_complexes for all
    op_notmany <- minres %>%
        decomplexify(columns = columns) %>%
        filter(!if_any(!ends_with("complex"), function(x) x %in% entity_2many)) %>%
        # recomplexify
        select(-columns) %>%
        rename_with(~gsub("_complex", "", .x), ends_with("complex")) %>%
        distinct()

    # homologous omnipath resource (with 1to1 alone)
    or_notmany <- op_notmany %>%
        # join back missing cols
        left_join(op_resource, by = columns) %>%
        .handle_complexes(symbols_dict = symbols_dict,
                          .missing_fun = .missing_fun,
                          columns=columns)

    ### Get interactions with many-many
    op_1many <- anti_join(minres,
                          op_notmany,
                          by = columns)

    # Recursively translate the 1many resource
    op_one2_many <- op_1many %>%
        liana::decomplexify(columns = columns) %>%
        group_by(across(ends_with("complex"))) %>%
        group_split()

    # On all one2many !!!
    suppressWarnings(pb <- dplyr::progress_estimated(length(op_one2_many)))
    or_many <- op_one2_many %>%
        map(function(op_row.decomp){
            if(verbose) pb$tick()$print()

            # create dictonary by row
            dicts_row <- .create_row_dict(symbols_dict,
                                          op_row.decomp,
                                          entity_2many,
                                          columns)

            # The tibble to be translated to all homologs
            op_row <- inner_join(op_row.decomp,
                                 op_1many,
                                 by = columns) %>%
                # recomplexify
                select(-columns) %>%
                rename_with(~gsub("_complex", "", .x),
                            ends_with("complex")) %>%
                distinct() %>%
                # join back remainder of cols
                left_join(op_resource,
                          by = columns)

            # Return all matches for all homologs
            or_row <- map(dicts_row, function(d){
                .handle_complexes(op_row,
                                  symbols_dict = d,
                                  .missing_fun = .missing_fun,
                                  columns = columns)
            }) %>%
                bind_rows() %>%
                distinct()
        }) %>%
        bind_rows()

    # Bind 1to1 and 1tomany
    or_resource <- bind_rows(or_notmany, or_many) %>%
        select(-ends_with("complex")) %>%
        select(!!columns, everything()) %>%
        distinct_at(c("source_genesymbol", "target_genesymbol"),
                    .keep_all = TRUE)

    return(or_resource)

}


#' Function to translate Human complexes to organism X homologs
#'
#' @param op_resource resource in the format of OmniPath/LIANA
#'
#' @param symbols_dict dictionary (named list) with human genesymbols and ortholog
#'  in a second species
#'
#' @param columns columns relevant for homology conversion
#'
#' @details Complexes cannot be joined directly, thus this function will
#' translate each subunit one by one.
#'
#' @noRd
.handle_complexes <- function(op_resource,
                              symbols_dict,
                              columns=columns,
                              .missing_fun = NULL){
    # decomplexify
    op_resource_decomplex <- op_resource %>%
        decomplexify(columns=columns)

    # translate subunits
    translated_subunits <- op_resource_decomplex %>%
        mutate(across(all_of(columns),
                      ~recode.character2(.x, !!!symbols_dict,
                                         .missing_fun = .missing_fun)))

    # Generate Dictionaries for complexes
    complex_dict <- map(columns,
                        ~.generate_complex_dict(translated_subunits, col=.x))


    # Bind all dictionaries
    dict <- pmap( # append multiple lists
        list(
            list(
                symbols_dict,
                flatten(complex_dict)
            )
        ), c) %>%
        flatten %>%
        flatten()

    # get orthologous resource
    op_ortholog <- op_resource %>%
        mutate(across(columns,
                      ~recode.character2(.x,
                                         !!!dict,
                                         .missing_fun = .missing_fun)))

    return(op_ortholog)
}


#' Modified `dplyr::recode` function
#'
#' @inheritParams dplyr::recode
#'
#' @param .missing_fun Function to modify any missing homologs/strings
#'  `NULL` by default and any missing values will be discarded.
#'  For example, one could be set it to `str_to_title` to format all symbols, or
#'  any other format in a scenario where a homolog dictionary is not available
#'  for the organism of interest.
#'
#' @details enables to modify unmatched genesymbols
#' @import stringr
#'
#' @keywords internal
recode.character2 <- function(.x,
                              ...,
                              .default = NULL,
                              .missing = NULL,
                              .missing_fun) {
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
        out <- dplyr:::replace_with(out,
                                    .x == nm,
                                    values[[nm]],
                                    paste0("`", nm, "`"))
        replaced[.x == nm] <- TRUE
    }

    .default <- dplyr:::validate_recode_default(.default, .x, out, replaced)

    if(!is.null(.missing_fun)){
        out <- dplyr:::replace_with(out,
                                    !replaced & !is.na(.x),
                                    exec(.missing_fun, .default),
                                    "`.default`")
        }

    out <- dplyr:::replace_with(out,
                                is.na(.x),
                                .missing,
                                "`.missing`")
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
.generate_complex_dict <- function(translated_subunits, col){
    col_complex <- str_glue("{col}_complex")

    filter(translated_subunits, str_detect(.data[[col_complex]], "_")) %>%
        dplyr::select(.data[[col_complex]], .data[[col]]) %>%
        distinct() %>%
        group_by(.data[[col_complex]]) %>%
        group_nest(keep = FALSE, .key = 'subunits') %>%
        mutate(translated_complex = map(subunits, function(sub){
            glue::glue_collapse(sub[[col]], sep = "_") %>%
                as.character()
        })) %>%
        dplyr::select(.data[[col_complex]], translated_complex) %>%
        unnest(translated_complex) %>%
        deframe()

}


#' Helper function to bind dictionaries
#'
#' @param main_entity entity with one-to-many mapping (ligand or receptor)
#' @param secondary_entity other entity (ligand or receptor)
#'
#' @details split by interaction, deframe and bind the secondary entity.
#' For example, if a ligand (`main_entity`) has one-to-many mapping to
#' multiple homologs, then main_entity will contain all homologs that match to
#' the ligand's subunits (or protein if not a complex). To this, we also attach
#' the genes of the receptor (`secondary_entity`) so that `generate_orthologs`
#' can translate both the ligand and receptor, and vice versa.
#'
#' @noRd
.bind_dicts <- function(main_entity, secondary_entity){
    main_entity %>%
        group_by_all() %>%
        group_split() %>%
        map(~deframe(bind_rows(.x, secondary_entity)))
}


#' Function to create a dictionary for one2many maps
#'
#' @param symbols_dict the human to mouse dict (should be a vect)
#' @param op_row.decomp decomp omni (we need all genes)
#' @param entity_2many entities that match to many homologs
#' @param columns columns
#'
#' @return a dictionary for the provided interaction
#'
#' @noRd
.create_row_dict <- function(symbols_dict,
                             op_row.decomp,
                             entity_2many,
                             columns){

    logical_row <- map_lgl(columns, function(col){
        any(op_row.decomp[[col]] %in% entity_2many)
    })

    dicts_row <- map(columns, function(col){
        symbols_dict[names(symbols_dict) %in% op_row.decomp[[col]]] %>%
            enframe(name = "genesymbol_source", value = "genesymbol_target")
    })

    # if only 1 col -> return dict
    if(length(columns)==1){
        symbols_dict[names(symbols_dict) %in% op_row.decomp[[columns]]] %>%
            enframe(name = "genesymbol_source", value = "genesymbol_target") %>%
            group_by_all() %>%
            group_split() %>%
            map(~deframe(.x))

    } else if(all(logical_row)){ # if all columns have 1-to-many
        c(.bind_dicts(dicts_row[[1]], dicts_row[[2]]),
          .bind_dicts(dicts_row[[2]], dicts_row[[1]]))
    } else{ # main entity is the one with 1-to-many
        .bind_dicts(main_entity = keep(dicts_row, logical_row) %>% pluck(1),
                    secondary_entity = discard(dicts_row, logical_row) %>% pluck(1))
    }
}


#' Helper function to get homologene dictionary
#'
#' @param entities Character vector: symbols of genes to be converted - this
#'     function returns a dictionary restricted to these genes.
#' @param target_organism Character or numeric: name or NCBI Taxonomy ID of the
#'     target organism.
#' @param mappings Character vector: control ambiguous mappings.
#'
#' @keywords internal
#'
#' @importFrom OmnipathR oma_pairwise
homology_dict <- function(
        entities,
        target_organism,
        id_type = "genesymbol",
        mappings = c("1:1", "1:m", "n:1", "n:m")
    ){

    oma_pairwise(
        organism_b = target_organism,
        id_type = id_type,
        mappings = mappings
    ) %>%
    # Limit to the universe of the resource
    filter(.data[["id_organism_a"]] %in% entities) %>%
    deframe

}


#' Deprecated call to generate_homologs
#'
#' @inheritDotParams generate_orthologs
#'
#' @export
generate_orthologs <- function(...){
    warning("`generate_orthologs` is deprecated, use `generate_homologs` instead")

    generate_homologs(...)
}

#' Helper function to show available organisms via OmnipathR's homologene resource
#'
#' @importFrom OmnipathR homologene_raw
#'
#' @export
show_homologene <- function(){
    homologene_raw() %>%
        pull(ncbi_taxid) %>%
        unique %>%
        tibble(ncbi_taxid = .,
               name = common_name(.),
               latin = latin_name(.)) %>%
        na.omit() %>%
        arrange(name) %>%
        print(n=nrow(.))
}

