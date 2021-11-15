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
#' @param omni_variants bool whether to get different OmniPath variants (e.g.
#' _full, based on ligrec resource quality quartile, or if only lig_rec)
#' @param lr_pipeline bool whether to format for lr_pipeline and remove
#'
#' duplicate LRs (mainly from composite OmniDB due to category (adhesion vs lr))
#' @return A list of OmniPath resources formatted according to the method pipes
#'
#' @importFrom magrittr %>%
#' @importFrom purrr pluck map
#' @importFrom rlang !!! exec
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
        filter(!(source %in% c("O75462", "Q13261"))) %>%
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
        distinct()


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
        select(-starts_with("interaction"))


    # CellChatDB Fix (append missing) ----
    ligrec$CellPhoneDB$interactions %<>%
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
reform_omni <- function(ligrec){
    map(ligrec, function(x) x %>%
            pluck("interactions") %>%
            distinct_at(.vars = c("source_genesymbol", # remove duplicate LRs
                                  "target_genesymbol"),
                        .keep_all = TRUE)) %>%
        append(list("Default" = NULL),
               .)
}



#' Retrieves intercellular interactions from OmniPath
#' @inheritParams omnipath_partners
#' @return A tibble with Intercell interactions from OmniPath
#'
#' @importFrom magrittr %>%
omnipath_intercell <- function(...){
    OmnipathR::import_intercell_network() %>%
        OmnipathR::filter_intercell_network(...)
}




#' Retrieves the interactions from one ligand-receptor resource
#' @inheritDotParams OmnipathR::import_post_translational_interactions
#' @inheritParams get_partners
#' @import tibble
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
omnipath_partners <- function(side, ...){

    causality <- list(ligand = 'trans', receptor = 'rec')

    OmnipathR::import_omnipath_intercell(
        causality = causality[[side]],
        scope = 'generic',
        source = 'composite',
        ...
    )

}


#' Shuffle OmniPath Intercell Resource
#' @param op_resource Intercell Resource to shuffled
#' @param .seed Value for set.seed
#'
#' @return A shuffled omnipath-formatted resource
#'
#' @import tibble
#' @noRd
shuffle_omnipath <- function(op_resource,
                             .seed = 1004){

    # library(BiRewire)
    set.seed(.seed)

    # make a vector proportional to the number of consensus directions
    stimul_num <- round(mean(op_resource$consensus_stimulation)/
                            mean(op_resource$consensus_direction) * 100000)
    directed_vector <- append(rep(1, stimul_num), rep(-1,100000 - stimul_num))

    op_prep <- op_resource %>%
        select(source_genesymbol, target_genesymbol,
               is_directed, is_stimulation, is_inhibition,
               consensus_direction, consensus_stimulation,
               consensus_inhibition) %>%
        mutate(mor = if_else(consensus_stimulation==1,
                             true=1,
                             if_else(consensus_inhibition==1, -1, 0))) %>%
        mutate(mor = if_else(mor==0, sample(directed_vector, 1), mor)) %>%
        select(source_genesymbol, mor, target_genesymbol) %>%
        distinct()

    # Induced bipartite and SIF builder
    op_dsg = BiRewire::birewire.induced.bipartite(
        op_prep,
        delimitators = list(negative = '-1',
                            positive = '1')
        )
    op_sif = BiRewire::birewire.build.dsg(
        op_dsg,
        delimitators = list(negative = '-1',
                            positive = '1')
        )

    # Rewire dsg
    random_dsg = BiRewire::birewire.rewire.dsg(dsg = op_dsg)
    random_sif = BiRewire::birewire.build.dsg(
        random_dsg,
        delimitators = list(negative = '-1',
                            positive = '1')
        )

    # Jacard dsg
    message(str_glue("Jaccard index between random and original resource: ",
                     {BiRewire::birewire.similarity.dsg(op_dsg, random_dsg)}))

    # format to OmniPath
    op_random <- random_sif %>%
        select(
            source,
            target,
            sign) %>%
        mutate(
            source_genesymbol = source,
            target_genesymbol = target,
            is_directed = 1,
            is_stimulation = if_else(sign==1, 1, 0),
            consensus_stimulation = if_else(sign==1, 1, 0),
            is_inhibition = if_else(sign==-1, 1, 0),
            consensus_inhibition = if_else(sign==-1, 1, 0),
            category_intercell_source = "ligand",
            category_intercell_target = "receptor",
            genesymbol_intercell_source = source_genesymbol,
            genesymbol_intercell_target = target_genesymbol,
            entity_type_intercell_target = "protein",
            sources = "RANDOM",
            references = "BiRewire",
            entity_type_intercell_source = "protein",
            entity_type_intercell_target
        ) %>%
        as_tibble()
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
                       "dip_url",    "sources","references", "curation_effort",
                       "n_references", "n_resources",
                       "category_intercell_source", "category_intercell_target")
            else .
        }

}



#' Function to Obtain the CellCall database
#'
#' @returns cellcall db converted to LIANA/OP format
get_cellcall <- function(){
    read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs.txt"), header = TRUE) %>%
        mutate(across(everything(), ~as.character(.x))) %>%
        # bind_rows(read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs_extended.txt"), header = TRUE) %>%
        #               mutate(across(everything(), ~as.character(.x)))) %>%
        select(source_genesymbol = Ligand_Symbol,
               target_genesymbol = Receptor_Symbol) %>%
        filter(source_genesymbol != target_genesymbol) %>%
        mutate(across(ends_with("genesymbol"), ~gsub("\\,", "\\_", .x))) %>%
        distinct() %>%
        filter(source_genesymbol != target_genesymbol) %>%
        # translate to Uniprot
        rowwise() %>%
        mutate(target = genesymbol_to_uniprot(target_genesymbol)) %>%
        mutate(source = genesymbol_to_uniprot(source_genesymbol)) %>%
        ungroup() %>%
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
        )
}


#' Helper Function to translate to UniProt
#' @param st any genesymbol string - to be separate by `_`
genesymbol_to_uniprot <- function(st){
    st.split <- as.vector(str_split(st, pattern = "_"))[[1]]
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys=st.split,
                          keytype = "SYMBOL",
                          columns = c("SYMBOL", "UNIPROT")) %>%
        arrange(desc(UNIPROT)) %>%
        distinct_at(c("SYMBOL"), .keep_all = TRUE) %>%
        pluck("UNIPROT") %>%
        glue::collapse(sep = "_")
}


#' Helper Function to get Missing Interactions from OG CellChatDB
get_cellchat_missing <- function(){
    CellChat::CellChatDB.human %>%
        pluck("interaction") %>%
        select(source_genesymbol = ligand,
               target_genesymbol = receptor) %>%
        mutate(across(everything(), ~str_to_upper(.x))) %>%
        mutate(across(everything(), ~gsub("*\\s..*" ,"" , .x))) %>%
        mutate(across(everything(), ~gsub("\\:" ,"_" , .x))) %>%
        # Obtain only ITGA1_ITGB1-interactions (they are missing from CellChatDB in Omni)
        filter(str_detect(target_genesymbol, "ITGA1_ITGB1")) %>%
        # translate to Uniprot
        rowwise() %>%
        mutate(target = genesymbol_to_uniprot(target_genesymbol)) %>%
        mutate(source = genesymbol_to_uniprot(source_genesymbol))
}

#' Helper function to obtain distinct transmitter and receiver lists
#' used in the resource comparison
#'
#' @param ligrec_list e.g. ligrec$OmniPath
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



#' Helper function to obtain unformatted CellCall DB from its source
#' @returns cellcall db formatted to OmniPath/LIANA
get_cellcall <- function(){
    read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs.txt"), header = TRUE) %>%
        mutate(across(everything(), ~as.character(.x))) %>%
        # we can also get extended interactions
        # bind_rows(read.delim(url("https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/new_ligand_receptor_TFs_extended.txt"), header = TRUE) %>%
        #               mutate(across(everything(), ~as.character(.x)))) %>%
        select(source_genesymbol = Ligand_Symbol,
               target_genesymbol = Receptor_Symbol) %>%
        filter(source_genesymbol != target_genesymbol) %>%
        mutate(across(ends_with("genesymbol"), ~gsub("\\,", "\\_", .x))) %>%
        distinct() %>%
        filter(source_genesymbol != target_genesymbol) %>%
        # translate to Uniprot
        rowwise() %>%
        mutate(target = genesymbol_to_uniprot(target_genesymbol)) %>%
        mutate(source = genesymbol_to_uniprot(source_genesymbol)) %>%
        ungroup() %>%
        # fill placeholders
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
        )
}
