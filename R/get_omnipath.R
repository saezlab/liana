#' Helper function that returns the name of each intercell resource in OmniPath
#'
#' @return A list of strings for each intercell resource in OmniPath
#'
#' @export
get_lr_resources <- function(){
    return(
        list(
            'Baccin2019',
            'CellCall',
            'CellChatDB',
            'Cellinker',
            'CellPhoneDB',
            'CellTalkDB',
            'connectomeDB2020',
            'EMBRACE',
            'Guide2Pharma',
            'HPMR',
            'ICELLNET',
            'iTALK',
            'Kirouac2010',
            'LRdb',
            'Ramilowski2015'
        )
    )
}

# only the ones different from the current defaults:
op_ic_quality_param <- list( # used for nodes
    resource = 'OmniPath', # this is just necessary in all the calls
    loc_consensus_percentile = 51,
    consensus_percentile = NULL
)

# used for interactions
op_ia_quality_param <- list(
    transmitter_topology = c('secreted',
                             'plasma_membrane_transmembrane',
                             'plasma_membrane_peripheral'),
    receiver_topology = c('plasma_membrane_transmembrane',
                          'plasma_membrane_peripheral'),
    min_curation_effort = 1,
    ligrecextra = FALSE
)


#' Function to get unfiltered intercell resources
#' For each resource and OmniPath variant compiles tables of ligands,
#' receptors and interactions
#'
#' @details calls on omnipath_intercell, intercell_connections, get_partners,
#' and intercell_connections
#'
#' @param lr_pipeline bool whether to format for lr_pipeline and remove
#'
#' duplicate LRs (mainly from composite OmniDB due to category (adhesion vs lr))
#' @return A list of OmniPath resources formatted according to the method pipes
#'
#' @importFrom magrittr %>%
#' @importFrom purrr pluck map
#' @importFrom rlang !!! exec
#'
#' @keywords internal
#'
#' @export
compile_ligrec <- function(lr_pipeline = TRUE){

    ligrec <-
        get_lr_resources() %>%
        map(function(resource){
            list(transmitters = get_ligands(resource),
                 receivers = get_receptors(resource),
                 interactions = intercell_connections(resource))
        }) %>%
        setNames(get_lr_resources()) %>%
        c(list(
            OmniPath = list(
                transmitters = exec(get_ligands, !!!op_ic_quality_param),
                receivers = exec(get_receptors, !!!op_ic_quality_param),
                interactions = exec(
                    intercell_connections,
                    !!!op_ic_quality_param,
                    !!!op_ia_quality_param
                )
            )
        ))

    # Format OmniPath ----
    ligrec$OmniPath$interactions %<>%
        filter(!(entity_type_intercell_source == "complex" |
                     entity_type_intercell_target == "complex"))  %>%
        # filter any mediators
        filter(!(str_detect(category_intercell_source, "cofactor")) &
                   !(str_detect(category_intercell_target, "cofactor")) &
                   !(str_detect(category_intercell_source, "ligand_regulator")))%>%
        # remove ambiguous/non-membrane associated receptor-receptor interactions
        # as well as others which seem to be misannotated (manually)
        filter(!(source %in% c("O75462", "Q13261", "P00533", "O00220",
                               "P06213", "P08254", "Q99835",
                               "Q9ULT6", "P06213", "Q13467", "P09619")
                 )) %>%
        # Filter KEA if it's the only curation
        filter(!(str_detect(sources, "KEA") & curation_effort==1)) %>%
        # filter any ion_channel/adp-associated interactions
        filter(parent_intercell_target != "ion_channel") %>%
        # interactions need to be reversed
        mutate(target_genesymbol_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                                       "CLEC2A", "CLEC2B", "CLEC2D"),
                                              source_genesymbol,
                                              target_genesymbol),
               source_genesymbol_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                                       "CLEC2A", "CLEC2B", "CLEC2D"),
                                              target_genesymbol,
                                              source_genesymbol)) %>%
        mutate(target_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                            "CLEC2A", "CLEC2B", "CLEC2D"),
                                              source,
                                              target),
               source_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                            "CLEC2A", "CLEC2B", "CLEC2D"),
                                   target,
                                   source)) %>%
        mutate(target = target_new,
               target_genesymbol = target_genesymbol_new,
               source = source_new,
               source_genesymbol = source_genesymbol_new) %>%
        dplyr::select(-ends_with("new")) %>%
        distinct() %>%
        select(-starts_with("plasma_membrane")) %>%
        select(source_genesymbol, target_genesymbol,
               source, target, everything())


    # Format CPDB ----
    ligrec$CellPhoneDB$interactions %<>%
        # check if any ambigous interactions (wrongly annotated ligands/receptors) exist
        rowwise() %>%
        unite(source, target, col = "interaction", remove = FALSE) %>%
        unite(target, source, col = "interaction2", remove = FALSE) %>%
        # identify duplicates
        mutate(dups = if_else(interaction %in% interaction2 |
                                  interaction2 %in% interaction,
                              TRUE,
                              FALSE)) %>%
        # ligands which are targets in OmniPath
        mutate(wrong_transitters = (source %in% ligrec$OmniPath$interactions$target)) %>%
        # receptors which are ligands in OmniPath
        mutate(wrong_receivers = (target %in% ligrec$OmniPath$interactions$source)) %>%
        # filter duplicates which are wrongly annotated
        filter(!(wrong_transitters & dups)) %>%
        filter(!(wrong_receivers & dups)) %>%
        # mismatched transmitters
        filter(!(target %in% c("P09917"))) %>%
        select(-starts_with("interaction"), -starts_with("wrong"), -dups)


    # CellChatDB Fix (append missing) ----
    ligrec$CellChatDB$interactions %<>%
        # append missing OG CellChatDB interactions
        bind_rows(get_cellchat_missing()) %>%
        mutate(across(where(is_double), ~replace_na(.x, 1))) %>%
        mutate(across(where(is_integer), ~replace_na(.x, 1))) %>%
        mutate(across(where(is_character), ~replace_na(.x, "placeholder")))


    # Obtain CellCall from source
    ligrec$CellCall$interactions <- get_cellcall()

    # Format to pipeline or not
    ligrec %<>% { if(lr_pipeline) reform_omni(.) else ligrec %>%
            map(function(ligrec_interactions){
                ## NOTE!!!
                # Obtain Transmitters and Receivers from the interactions
                ligrec_interactions %<>%
                    assign_ligrecs()
                })
        }

    return(ligrec)
}



#' Helper Function to Reformat ligrec for LR Pipeline
#' @param ligrec OmniPath list returned by compile_ligrec
#' @return A list of OmniPath resources, including OmniPath composite DB,
#' A reshuffled OmniPath, and a Default with NULL ( tool pipelines run
#' using their default resource)
#' @importFrom purrr pluck map
#' @importFrom dplyr distinct_at
#'
#' @noRd
reform_omni <- function(ligrec){
    map(ligrec, function(x) x %>%
            pluck("interactions") %>%
            distinct_at(.vars = c("source_genesymbol", # remove duplicate LRs
                                  "target_genesymbol"),
                        .keep_all = TRUE)) %>%
        append(list("Default" = NULL,
                    "Consensus" = get_curated_omni()),
               .)
}



#' Retrieves intercellular interactions from OmniPath
#' @inheritParams omnipath_partners
#' @inheritDotParams OmnipathR::filter_intercell_network
#' @return A tibble with Intercell interactions from OmniPath
#'
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_intercell <- function(...){
    OmnipathR::import_intercell_network() %>%
        OmnipathR::filter_intercell_network(...)
}




#' Retrieves the interactions from one ligand-receptor resource
#'
#' @inheritDotParams OmnipathR::import_post_translational_interactions
#' @inheritParams get_partners
#' @import tibble
#'
#' @noRd
intercell_connections <- function(resource, ...){

    if(resource == 'OmniPath'){

        return(omnipath_intercell(...))

    }

    OmnipathR::import_post_translational_interactions(
        resource = resource,
        ...
    ) %>%
        as_tibble() %>%
        mutate(category_intercell_source = "ligand",
               category_intercell_target = "receptor")

}


#' Retrieves ligands from one ligand receptor resource
#' @inheritDotParams intercell_connections
#' @inheritParams get_partners
#'
#' @noRd
get_ligands <- function(resource, ...){

    get_partners(side = 'ligand', resource = resource, ...)

}


#' Retrieves receptors from one ligand-receptor resource
#'
#' @inheritDotParams intercell_connections
#'
#' @noRd
get_receptors <- function(resource, ...){

    get_partners(side = 'receptor', resource = resource, ...)

}



#' Retrieves intercellular communication partners (ligands or receptors) from
#' one ligand-receptor resource.
#' @inheritParams omnipath_partners
#' @param resource Name of current resource (taken from get_lr_resources)
#' @param ... Inherit dot params from \link{OmnipathR::omnipath_intercell}
#' @importFrom rlang sym !!!
#' @importFrom magrittr %>%
get_partners <- function(side, resource, ...){

    if(resource == 'OmniPath'){
        return(omnipath_partners(side = side, ...))
    }

    id_cols <- `if`(
        side == 'ligand',
        syms(c('source', 'source_genesymbol')),
        syms(c('target', 'target_genesymbol'))
    )
    up_col <- `if`(
        side == 'ligand',
        sym('source'),
        sym('target')
    )
    gs_col <- `if`(
        side == 'ligand',
        sym('source_genesymbol'),
        sym('target_genesymbol')
    )

    intercell_connections(resource, ...) %>%
        select(!!!id_cols) %>%
        distinct() %>%
        rename(uniprot = !!up_col, genesymbol = !!gs_col)

}




#' Retrieves intercellular communication partners (transmitters or receivers)
#' from OmniPath
#'
#' @param side 'ligand' (trans), 'receptor' (rec) or 'both' (both short or
#'     long notation can be used)
#' @inheritDotParams OmnipathR::import_omnipath_intercell
#'
#' @importFrom OmnipathR import_omnipath_intercell
#' @importFrom magrittr %>%
#'
#' @keywords internal
omnipath_partners <- function(side, ...){

    causality <- list(ligand = 'trans', receptor = 'rec')

    OmnipathR::import_omnipath_intercell(
        causality = causality[[side]],
        scope = 'generic',
        source = 'composite',
        ...
    )

}


#' Function to Generate OmniPath versions
#'
#' @param remove_complexes whether to remove complexes
#' @param simplify whether to simplify according to the mandatory columns needed by different methods in `liana`
#' @inheritDotParams OmnipathR::filter_intercell_network
#'
#' @export
#'
#' @return An OmniPath resource
#'
#' @noRd
generate_omni <- function(remove_complexes=TRUE,
                          simplify = TRUE,
                          ...){
    OmnipathR::import_intercell_network() %>%
        {
            if(remove_complexes)
                filter(., !(entity_type_intercell_source == "complex" |
                                entity_type_intercell_target == "complex"))
            else .

        } %>%
        OmnipathR::filter_intercell_network(
            simplify = FALSE,
            ...
        )  %>%
        distinct_at(.vars = c("source_genesymbol", # remove duplicate LRs
                              "target_genesymbol"),
                    .keep_all = TRUE) %>%
        {
            if(simplify)
                select(.,
                       "source", "target", "source_genesymbol", "target_genesymbol",
                       "is_directed", "is_stimulation", "is_inhibition",
                       "consensus_direction","consensus_stimulation", "consensus_inhibition",
                       "sources","references", "curation_effort",
                       "n_references", "n_resources",
                       "category_intercell_source", "category_intercell_target")
            else .
        }

}



#' Function to Obtain the CellCall database
#'
#' @keywords internal
#'
#' @returns cellcall db converted to LIANA/OP format
get_cellcall <- function(){
    # Get UniProt Query DB
    cellcall <- read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs.txt"), header = TRUE) %>%
        mutate(across(everything(), ~as.character(.x))) %>%
        # we can also get extended interactions
        # bind_rows(read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs_extended.txt"), header = TRUE) %>%
        #               mutate(across(everything(), ~as.character(.x)))) %>%
        dplyr::select(Ligand_ID,
                      Receptor_ID) %>%
        filter(Ligand_ID != Receptor_ID) %>%
        mutate(across(everything(), ~gsub("\\,", "\\_", .x))) %>%
        distinct()

    up <- UniProt.ws::UniProt.ws(taxId=9606)
    # Get Dict
    up_dict <- get_up_dict(cellcall, up)

    # convert to lists for recode
    up_dict_uniprot <- up_dict %>%
        dplyr::select(GENEID, uniprot) %>%
        mutate(across(everything(), ~gsub("\\s", "", .x))) %>%
        deframe() %>%
        as.list()
    up_dict_symbol <- up_dict %>%
        dplyr::select(GENEID, genesymbol) %>%
        mutate(genesymbol = gsub("*\\s..*" ,"" , genesymbol)) %>%
        mutate(across(everything(), ~gsub("\\s", "", .x))) %>%
        deframe() %>%
        as.list()

    # translate to Uniprot and format
    cellcall %>%
        rowwise() %>%
        mutate(source = geneid_to_uniprot(Ligand_ID, up_dict_uniprot)) %>%
        mutate(target = geneid_to_uniprot(Receptor_ID, up_dict_uniprot)) %>%
        mutate(source_genesymbol = geneid_to_uniprot(Ligand_ID, up_dict_symbol)) %>%
        mutate(target_genesymbol = geneid_to_uniprot(Receptor_ID, up_dict_symbol)) %>%
        ungroup() %>%
        select(-c(Ligand_ID, Receptor_ID)) %>%
        mutate(is_directed = 1,
               is_stimulation = 1,
               is_inhibition = 1,
               consensus_direction = 1,
               consensus_stimulation = 1,
               consensus_inhibition = 1,
               dip_url = "placeholder",
               sources = "placeholder",
               references = "placeholder",
               curation_effort = "placeholder",
               n_references = 1,
               n_resources = 1,
               category_intercell_source = "placeholder",
               category_intercell_target = "placeholder"
        ) %>%
        mutate(across(c("source", "target",
                        "source_genesymbol",
                        "target_genesymbol"),
                      ~as.character(.x)))
}


#' Helper Function to translate to UniProt
#' @param st any genesymbol string - to be separate by `_`
#'
#' @keywords internal
genesymbol_to_uniprot <- function(st){
    st.split <- as.vector(str_split(st, pattern = "_"))[[1]]
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys=st.split,
                          keytype = "SYMBOL",
                          columns = c("SYMBOL", "UNIPROT", "EVIDENCE")) %>%
        arrange(desc(UNIPROT)) %>%
        distinct_at(c("SYMBOL"), .keep_all = TRUE) %>%
        pluck("UNIPROT") %>%
        glue::collapse(sep = "_")
}


#' Helper function to get UniProt dictionary
#'
#' @param ligrec_res ligand_receptor resource to translate
#' @param up uniprot db to be queried
#' @param key_column1 name of the ligand  column
#' @param key_column2 name of the receptor column
#'
#' @keywords internal
get_up_dict <- function(ligrec_res,
                        up,
                        key_column1 = "Ligand_ID",
                        key_column2 = "Receptor_ID"){

    keys <- unlist(union(str_split(ligrec_res[[key_column1]], pattern = "_"),
                         str_split(ligrec_res[[key_column2]], pattern = "_")))

    up_dict <- UniProt.ws::select(up, keytype = c("GENEID"),
                                  columns = c("UNIPROTKB", "GENES","REVIEWED"),
                                  keys = keys) %>%
        filter(REVIEWED == "reviewed") %>%
        dplyr::select(GENEID,
                      uniprot = UNIPROTKB,
                      genesymbol = GENES) %>%
        # Keep only first gene symbol (i.e. official one)
        tibble() %>%
        mutate(GENEID = gsub("\\s", "", GENEID))

    return(up_dict)
}

#' Helper Function to translate to UniProt
#'
#' @param st any genesymbol string - to be separate by `_`
#' @param dict dictionary with genesymbols and uniprot IDs
#'
#' @noRd
geneid_to_uniprot <- function(st,
                              dict){
    st.split <- str_split(st, pattern = "_")[[1]]

    tryCatch(
        {
            map(st.split, function(spl){
                recode(as.character(gsub("\\s", "", spl)), !!!dict)
            }) %>%
                glue::glue_collapse(sep = "_")

        },
        error = function(cond){
            message(str_glue("{st} had no match!"))
            return(NA)
        }
    )
}



#' Helper Function to get Missing Interactions from OG CellChatDB
#'
#' @keywords internal
get_cellchat_missing <- function(){
    CellChat::CellChatDB.human %>%
        pluck("interaction") %>%
        select(source_genesymbol = ligand,
               target_genesymbol = receptor) %>%
        mutate(across(everything(), ~stringr::str_to_upper(.x))) %>%
        mutate(across(everything(), ~gsub("*\\s..*" ,"" , .x))) %>%
        mutate(across(everything(), ~gsub("\\:" ,"_" , .x))) %>%
        # Obtain only ITGA1_ITGB1-interactions (they are missing from CellChatDB in Omni)
        filter(str_detect(target_genesymbol, "ITGA1_ITGB1")) %>%
        # translate to Uniprot
        rowwise() %>%
        mutate(target = genesymbol_to_uniprot(target_genesymbol)) %>%
        mutate(source = genesymbol_to_uniprot(source_genesymbol)) %>%
        ungroup()
}

#' Helper function to obtain distinct transmitter and receiver lists
#' used in the resource comparison
#'
#' @param ligrec_list e.g. ligrec$OmniPath
#'
#' @keywords internal
assign_ligrecs <- function(ligrec_list){
    ligrec_list[["transmitters"]] <- ligrec_list$interactions %>%
        select(genesymbol = source_genesymbol,
               uniprot = source) %>%
        distinct()
    ligrec_list[["receivers"]] <- ligrec_list$interactions %>%
        select(genesymbol = target_genesymbol,
               uniprot = target) %>%
        distinct()

    return(ligrec_list)
}



#' Helper Function to check if there are dissociated entities which also exist as
#' complexes.
#'
#' @details We count the times that a ligand (check_entity) exists in a combination
#' with the same receptor (anchor_entity)
#'
#' @param complex_omni an OmniPath resource with complexes
#' @param check_entity the entity to be check for duplicates
#' @param anchor_entity the anchor entity with which we check for duplicates
#'
#' @keywords internal
check_if_dissociated <- function(complex_omni, check_entity, anchor_entity){
    check.complex <- str_glue("{check_entity}_complex")
    anchor.complex <- str_glue("{anchor_entity}_complex")

    complex_omni %>%
        liana::decomplexify() %>%
        select(-anchor_entity) %>%
        distinct() %>%
        group_by(.data[[check_entity]], .data[[anchor.complex]]) %>%
        # count the number that the entity exists with the same target
        mutate(counter = n()) %>%
        select(check_entity,
               check.complex, anchor.complex, counter) %>%
        arrange(desc(counter)) %>%
        # if exists more than once and is not a complex, then this
        # check_entity and anchor_entity is duplicated
        # hence we can remove any non-complex check_entity (since it already exists)
        filter(counter > 1) %>%
        filter(!str_detect(.data[[check.complex]], "_")) %>%
        ungroup() %>%
        # recomplexify
        select(-check_entity) %>%
        dplyr::rename({{check_entity}} := check.complex,
                      {{anchor_entity}} := anchor.complex)
}



#' Function to Generate the Curated (Default) LIANA resource
#'
#' @param curated_resources the curated resources from which we wish to obtain interactions.
#' By default, it includes interactions curated in the context of CCC from CellPhoneDB,
#' CellChat, ICELLNET, connectomeDB, CellTalkDB, and SignaLink.
#'
#' @details Here, we define the curated resources as those that are defined as manually
#' or expert curated in the context of cell-cell communication.
#' Albeit, "Guide2Pharma", "HPMR", and "Kirouac2010" are also such resources the remainder
#' of the resources used to generate Omnipath, use those as sources.
#' Hence, we assume that the second round of manual curation done in subsequent, more recently published
#' resources would already contain the high quality interactions of the aforementioned 3.
#' We also omit Cellinker, as it results in a large mount of ambigous interactions, but
#' one could consider adding it to the list of curated resources.
#'
#' @return a curated OmniPath resource formatted for LIANA
#'
#' @export
#'
get_curated_omni <- function(curated_resources = c("CellPhoneDB",
                                                   "CellChatDB",
                                                   "ICELLNET",
                                                   "connectomeDB2020",
                                                   "CellTalkDB")){


    # import the OmniPathR intercell network component
    ligrec <- OmnipathR::import_intercell_network(resources = curated_resources)

    # 1. Distinct and remove odd ligands
    complex_omni <- ligrec %>%
        filter(!category_intercell_source %in% c("activating_cofactor",
                                                 "ligand_antagonist")) %>%
        # we keep only the distinct LRs
        # (technically any combination of subunits can exist at this stage)
        distinct_at(c("source_genesymbol", "target_genesymbol"), .keep_all = TRUE) %>%
        # we then decomplexify (or split all complexes into subunits)
        liana::decomplexify(columns = c("source_genesymbol", "target_genesymbol"))

    # 2. Identify resulting duplicates from decomplexify
    # (keep the false/duplicate interactions alone)
    duplicated_lrs_only <- complex_omni %>%
        group_by(source_genesymbol, target_genesymbol) %>%
        # Iterative counter/ticker
        mutate(number = 1) %>%
        mutate(ticker = cumsum(number)) %>%
        filter(ticker > 1) %>%
        # which ones are complexes (they are unique)
        mutate(complex_flag = str_detect(source, pattern = "COMPLEX") |
                   str_detect(target, pattern = "COMPLEX")) %>%
        # remove any that are duplicated and not complexes
        filter(!complex_flag)

    # We anti join the false interactions
    complex_omni %<>%
        anti_join(duplicated_lrs_only)

    # 3. Remove duplicated complexes (introduced by expanding them via liana_decomplexify)
    complex_omni %<>%
        distinct_at(c("source_genesymbol_complex", "target_genesymbol_complex"), .keep_all = TRUE) %>%
        select(-c("source_genesymbol", "target_genesymbol")) %>%
        dplyr::rename(source_genesymbol = source_genesymbol_complex,
                      target_genesymbol = target_genesymbol_complex)

    # 4. Filter by localisation
    complex_omni_temp <- complex_omni %>%
        OmnipathR::filter_intercell_network(loc_consensus_percentile = 33,
                                            consensus_percentile = NULL,
                                            transmitter_topology = c('secreted',
                                                                     'plasma_membrane_transmembrane',
                                                                     'plasma_membrane_peripheral'),
                                            receiver_topology = c('plasma_membrane_transmembrane',
                                                                  'plasma_membrane_peripheral'),
                                            min_curation_effort = 0,
                                            min_resources = 1,
                                            min_references = 0,
                                            min_provenances = 0,
                                            simplify = TRUE) %>%
        select(-starts_with("is"))

    # since complexes are *penalized in terms of loc_consensus due to the fact
    # that they have multiple subunits the localisation annotations, hence
    # we don't filter them based on locallisation alone
    complex_omni_diff <- anti_join(complex_omni, complex_omni_temp) %>%
        select(colnames(complex_omni_temp)) %>%
        filter(!str_detect(source_genesymbol, "_") & !str_detect(target_genesymbol, "_"))

    # Filter out interactions (between non-proteins) whose localisation is
    # not appropriately assigned and/or are not referenced to a specific resource
    complex_omni %<>%
        select(colnames(complex_omni_temp)) %>%
        anti_join(complex_omni_diff)


    # 5. Get *curated interactions alone
    # Here we collected the curated interactions from all resources which contained any.
    # We defined curated interactions as those which contain a corresponding PubMed ID
    # Note that  we check for PubMed ID interactions solely as another level of certainty
    ligrec_curated <- OmnipathR::curated_ligand_receptor_interactions(
        curated_resources = curated_resources,
        cellphonedb = TRUE,
        cellinker = TRUE,
        talklr = FALSE,
        signalink = TRUE) %>%
        liana:::decomplexify(columns = c("source_genesymbol", "target_genesymbol"))

    # Check which curated interactions are missing from OmniPath
    curated_not_in_omni <- ligrec_curated %>%
        anti_join(complex_omni %>%
                      liana:::decomplexify(columns = c("source_genesymbol",
                                                       "target_genesymbol")),
                  by=c("source_genesymbol", "target_genesymbol")) %>%
        distinct_at(c("source_genesymbol", "target_genesymbol"), .keep_all = TRUE) %>%
        # Filter any ligands that are receptors in the localisation-filtered OmniPath
        filter(!(source_genesymbol %in% complex_omni$target_genesymbol)) %>%
        # and vice versa
        filter(!(target_genesymbol %in% complex_omni$source_genesymbol)) %>%
        # re-complexify
        select(-c("source_genesymbol", "target_genesymbol")) %>%
        dplyr::rename(source_genesymbol = source_genesymbol_complex,
                      target_genesymbol = target_genesymbol_complex) %>%
        # remove unneeded columns
        select(source, target, source_genesymbol, target_genesymbol, sources, references)

    # Bind the missing curated (Pubmed supported) and localisation-corrected interactions
    complex_omni %<>%
        bind_rows(curated_not_in_omni)

    # Remove left-over Redundant Subunits from Complex Ligands
    # (mostly the case when both the ligand and receptor are complexes)
    complex_omni %<>%
        # remove redundant ligands
        anti_join(check_if_dissociated(complex_omni,
                                       check_entity = "source_genesymbol",
                                       anchor_entity = "target_genesymbol"),
                  by = c("source_genesymbol", "target_genesymbol")) %>%
        # remove redundant receptors
        anti_join(check_if_dissociated(complex_omni,
                                       check_entity = "target_genesymbol",
                                       anchor_entity = "source_genesymbol"),
                  by = c("source_genesymbol", "target_genesymbol")) %>%
        mutate(across(everything(), ~replace_na(.x, "")))

    # 6. Check for duplicated/ambigous interactions
    # These include examples such as L1_R1 and R1_L1 ->
    # in general, we want to remove these,
    # unless they are in the context of cell-cell adhesion.
    duplicated_omni <- complex_omni %>%
        # check if any ambiguous/duplicated interactions exist
        rowwise() %>%
        unite(source, target, col = "interaction", remove = FALSE) %>%
        unite(target, source, col = "interaction2", remove = FALSE) %>%
        # identify duplicated interactions
        mutate(dups = if_else(interaction %in% interaction2 |
                                  interaction2 %in% interaction,
                              TRUE,
                              FALSE)) %>%
        # ligands which are targets in OmniPath
        mutate(ambigous_transmitters = (source %in% target)) %>%
        # filter to ambiguously annotated transmitters from duplicates
        # interactions alone
        filter((ambigous_transmitters & dups))

    wrong_transmitters <- c("ADGRE5", "CD160",
                            "CD226", "EGFR",
                            "TNFRSF18", "CTLA4",
                            "KLRB1", "KLRF1", "KLRF2",
                            "PTPRC", "PVR", "SIGLEC1",
                            "SIGLEC9", "TNFRSF14",
                            "ITGAD_ITGB2",
                            "ITGA4_ITGB1", "ITGA9_ITGB1",
                            "ITGA4_ITGB7",
                            "TYK2", "SYK",
                            "MT-RNR2",
                            "IL13_IL13RA1_IL4R",
                            "IL22_IL22RA1",
                            "IL18BP"

    )

    # # Identify wrongly annotated interactions
    duplicated_omni %<>%
        filter(source_genesymbol %in% wrong_transmitters)

    # Blocklist certain receivers
    block_receivers <- c("IFNG_IFNGR1", # include a ligand in the complex
                         "CNTN2_CNTNAP2",
                         "IL2_IL2RA_IL2RB_IL2RG",
                         "IL15_IL15RA_IL2RB_IL2RG",
                         "IL6_IL6R_IL6ST",
                         "IL1B_IL1R1_IL1RAP",
                         "IL1B_IL1R2_IL1RAP",
                         "IFNA2_IFNAR1_IFNAR2",
                         "ACVR1C_ACVR2B_CFC1",
                         "CSF2_CSF2RA_CSF2RB",
                         "GP1BA_GP1BB_GP5_GP9")


    # Remove those from our curated omnipath
    complex_omni %<>%
        anti_join(duplicated_omni, by=c("source_genesymbol",
                                        "target_genesymbol")) %>%
        filter(!target_genesymbol %in% block_receivers) %>%
        filter(!source_genesymbol %in% wrong_transmitters) %>%
        # # Recomplexify CD8 complex?
        # mutate(target_genesymbol = if_else(target_genesymbol=="CD8A" |
        #                                        target_genesymbol=="CD8B",
        #                                    "CD8A_CD8B",
        #                                    target_genesymbol)) %>%
        distinct_at(.vars=c("source_genesymbol", "target_genesymbol"),
                    .keep_all=TRUE)

    ## Recode odd aliases
    complex_omni %<>%
        mutate(source_genesymbol =
                   recode(source_genesymbol,
                          !!!list(
                              "C4B_2"="C4B"
                              )
                          )
               )

    # Return the final Curated OmniPath Resource
    return(complex_omni)
}

