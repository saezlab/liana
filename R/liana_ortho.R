
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
#'
#' @param .missing_fun approach to handle missing interactions. By default
#' set to `NULL` which would mean that any interactions without a homology
#' match will be filtered. This can be set to e.g. `str_to_title` when working
#' with murine symbols. Then if a gene has no matched homolog, instead of
#' discarding it, the `.missing_fun` will be used to format the name from human.
#' Hence, increasing the number of matches, but likely introducing some
#' mismatches.
#'
#' @param symbols_dict `NULL` by default, then `get_homologene_dict` is called
#' to generate a dictionary from OmniPathR's homologene resource. Alternatively,
#' one can pass their own symbols_dictionary.
#'
#' @param source name of the source (ligand) column
#'
#' @param target name of the target (receptor) column
#'
#' @param verbose logical for verbosity
#'
#' @return a converted ligand-receptor resource
#'
#' @export
generate_homologs <- function(op_resource,
                              target_organism,
                              max_homologs = 5,
                              .missing_fun = NULL,
                              symbols_dict = NULL,
                              source = "source_genesymbol",
                              target = "target_genesymbol",
                              verbose = TRUE){

    # Break down resource into subunit
    op_resource.decomp <- op_resource %>%
        liana::decomplexify(columns = c(source, target))

    # Get union of symbols
    entities <- union(op_resource.decomp[[source]],
                      op_resource.decomp[[target]])

    if(is.null(symbols_dict)){
        symbols_dict <- get_homologene_dict(entities,
                                            target_organism)
    }

    # Remove any missing antities
    if(is.null(.missing_fun)){

        # All missing entities
        missing_entities <- setdiff(union(op_resource.decomp$source_genesymbol,
                                          op_resource.decomp$target_genesymbol),
                                    names(symbols_dict))

        liana_message(
            str_glue("Entries without homologs:
                     {paste(missing_entities, collapse = '; ')}"),
            verbose = verbose
        )

        # Keep only interactions for which all proteins (incl. subunits) have a matching homologue
        op_resource.decomp %<>%
            mutate(source_present = !(source_genesymbol %in% missing_entities)) %>%
            mutate(target_present = !(target_genesymbol %in% missing_entities)) %>%
            mutate(lr_present = source_present * target_present) %>%
            mutate(source_genesymbol=source_genesymbol_complex,
                   target_genesymbol=target_genesymbol_complex) %>%
            group_by(source, target, source_genesymbol, target_genesymbol) %>%
            summarise(all_present = mean(lr_present), .groups = "keep") %>%
            # only keep those that are present
            filter(all_present < 1)

        # Remove interactions without matches
        op_resource <- anti_join(op_resource,
                                 op_resource.decomp,
                                 by = c("source", "target",
                                        "source_genesymbol", "target_genesymbol"))
    }

    # Obtain genes with multiple matches
    entity_2many <- symbols_dict %>%
        enframe(name = "genesymbol_source", value = "genesymbol_target") %>%
        group_by(genesymbol_source) %>%
        count(name = "n_match") %>%
        filter(n_match > 1 & n_match <= max_homologs) %>%
        pull(genesymbol_source)

    liana_message(
        str_glue("One-to-many homolog matches: {paste(entity_2many, collapse = '; ')}"),
        verbose = verbose
    )

    ### IF NOT many2many .handle_complexes for all
    op_notmany <- op_resource %>%
        liana::decomplexify() %>%
        filter(!(source_genesymbol %in% entity_2many),
               !(target_genesymbol %in% entity_2many))

    # homologous omnipath resource (with 1to1 alone)
    or_notmany <- op_notmany %>%
        # recomplexify
        mutate(source_genesymbol = source_genesymbol_complex,
               target_genesymbol = target_genesymbol_complex) %>%
        distinct() %>%
        .handle_complexes(symbols_dict = symbols_dict,
                          .missing_fun = .missing_fun)


    ### Get interactions with many-many
    op_1many <- anti_join(op_resource, or_notmany,
                          by = c("source", "target",
                                 "source_genesymbol", "target_genesymbol",
                                 "category_intercell_source",
                                 "database_intercell_source",
                                 "category_intercell_target",
                                 "database_intercell_target",
                                 "sources", "references"))

    # Recursively translate the 1many resource
    op_one2_many <- op_1many %>%
        liana::decomplexify() %>%
        group_split(source_genesymbol_complex,
                    target_genesymbol_complex)

    # On all one2many !!!
    suppressWarnings(pb <- dplyr::progress_estimated(length(op_one2_many)))
    or_many <- op_one2_many %>%
        map(function(op_row.decomp){
            if(verbose) pb$tick()$print()

            dicts_row <- .create_row_dict(symbols_dict,
                                          op_row.decomp,
                                          entity_2many)


            # The tibble to be translated to all homologs
            op_row <- inner_join(op_row.decomp, op_1many,
                                 by = c("source", "target",
                                        "category_intercell_source",
                                        "database_intercell_source",
                                        "category_intercell_target",
                                        "database_intercell_target",
                                        "sources", "references",
                                        "source_genesymbol", "target_genesymbol")) %>%
                mutate(source_genesymbol=source_genesymbol_complex,
                       target_genesymbol=target_genesymbol_complex) %>%
                distinct()

            # Return all matches for all homologs
            or_row <- map(dicts_row, function(d){
                .handle_complexes(op_row,
                                  symbols_dict = d,
                                  .missing_fun = .missing_fun)
            }) %>%
                bind_rows() %>%
                distinct()
        }) %>%
        bind_rows() %>%
        select(-ends_with("complex"))

    # Bind 1to1 and 1tomany
    or_resource <- bind_rows(or_notmany, or_many) %>%
        select(-ends_with("complex")) %>%
        select(!!source, !!target, everything())

    return(or_resource)

}

#' Function to translate Human complexes to organism X homologs
#'
#' @param op_resource resource in the format of OmniPath/LIANA
#'
#' @param symbols_dict dictionary (named list) with human genesymbols and ortholog
#'  in a second species
#'
#' @details Complexes cannot be joined directly, thus this function will
#' translate each subunit one by one.
#'
#' @noRd
.handle_complexes <- function(op_resource,
                              symbols_dict,
                              .missing_fun = NULL){
    # decomplexify
    op_resource_decomplex <- op_resource %>%
        decomplexify()

    # translate subunits
    translated_subunits <- op_resource_decomplex %>%
        mutate(across(ends_with("genesymbol"),
                      ~recode.character2(.x, !!!symbols_dict,
                                         .missing_fun = .missing_fun)))

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
        group_by(genesymbol_source, genesymbol_target) %>%
        group_split() %>%
        map(~deframe(bind_rows(.x, secondary_entity)))
}

#' Function to create a dictionary for one2many maps
#'
#' @param symbols_dict the human to mouse dict (should be a vect)
#' @param op_row.decomp decomp omni (we need all genes)
#' @param entity_2many entities that match to many homologs
#'
#' @return a dictionary for the provided interaction
#'
#' @noRd
.create_row_dict <- function(symbols_dict,
                             op_row.decomp,
                             entity_2many
                             ){
    # do the ligand or receptor posses homologs
    is.l.2many <- any(op_row.decomp$source_genesymbol %in% entity_2many)
    is.r.2many <- any(op_row.decomp$target_genesymbol %in% entity_2many)

    # Ligands
    dicts_row.l <-
        symbols_dict[names(symbols_dict) %in% op_row.decomp$source_genesymbol] %>%
        enframe(name = "genesymbol_source", value = "genesymbol_target")
    # Receptors
    dicts_row.r <-
        symbols_dict[names(symbols_dict) %in% op_row.decomp$target_genesymbol] %>%
        enframe(name = "genesymbol_source", value = "genesymbol_target")

    if(all(is.l.2many, is.r.2many)){
        c( # both
            .bind_dicts(dicts_row.l, dicts_row.r),
            .bind_dicts(dicts_row.r, dicts_row.l)
        )

    } else if(is.l.2many){
        .bind_dicts(dicts_row.l, dicts_row.r)
    } else if(is.r.2many){
        .bind_dicts(dicts_row.r, dicts_row.l)
    }
}



#' Helper function to get homologene dictionary
#'
#' @param entities genes to be converted
#'
#' @param target_organism target organism (obtain tax id from `show_homologene`)
#'
#' @keywords internal
#'
#' @importFrom OmnipathR homologene_download
get_homologene_dict <- function(entities,
                                target_organism){

    # Load homology geneset
    hg_gs <- homologene_download(target = !!target_organism,
                                 source = 9606L, # always human
                                 id_type = "genesymbol") %>%
        select(-hgroup) %>%
        # Limit to the universe of the resource
        filter(genesymbol_source %in% entities)

    # Convert to dictionary
    return(hg_gs %>% deframe())
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

