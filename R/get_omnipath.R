#' Helper function that returns the name of each intercell resource in OmniPath
#'
#' @return A list of strings for each intercell resource in OmniPath
#' @export
get_lr_resources <- function(){

    return(
        list(
            'CellChatDB',
            'CellPhoneDB',
            'Ramilowski2015',
            'Baccin2019',
            'LRdb',
            'Kirouac2010',
            'ICELLNET',
            'iTALK',
            'EMBRACE',
            'HPMR',
            'Guide2Pharma',
            'connectomeDB2020',
            'talklr',
            'CellTalkDB'
        )
    )

}

# only the ones different from the current defaults:
op_ic_quality_param <- list( # used for nodes
    resource = 'OmniPath', # this is just necessary in all the calls
    loc_consensus_percentile = 30,
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
#' @details calls on omnipath_intercell, intercell_connections, get_partners,
#' and intercell_connections
#' @param omni_variants bool whether to get different OmniPath variants (e.g.
#' _full, based on ligrec resource quality quartile, or if only lig_rec)
#' @param lr_pipeline bool whether to format for lr_pipeline and remove
#' duplicate LRs (mainly from composite OmniDB due to category (adhesion vs lr))
#' @return A list of OmniPath resources formatted according to the method pipes
#' @importFrom magrittr %>%
#' @importFrom purrr pluck map
#' @importFrom rlang !!! exec
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

    # Keep only nodes that are part of the interactions
    ligrec$OmniPath$receivers %<>%
        dplyr::filter(genesymbol %in% ligrec$OmniPath$interactions$target_genesymbol) %>%
        distinct_at(.vars="genesymbol", .keep_all = TRUE)
    ligrec$OmniPath$transmitters %<>%
        dplyr::filter(genesymbol %in% ligrec$OmniPath$interactions$source_genesymbol) %>%
        distinct_at(.vars="genesymbol", .keep_all = TRUE)

    ligrec %<>% {if(lr_pipeline) reform_omni(.) else .}

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
        append(list("Reshuffled" = shuffle_omnipath(.$connectomeDB2020),
                    "Default" = NULL),
               .)
}



#' Retrieves intercellular interactions from OmniPath
#' @inheritParams omnipath_partners
#' @return A tibble with Intercell interactions from OmniPath
#'
#' @importFrom magrittr %>%
# #'@importFrom OmnipathR import_intercell_network filter_intercell_network
#' @import OmnipathR
omnipath_intercell <- function(...){
    import_intercell_network(entity_types = 'protein') %>%
        filter_intercell_network(...)
}




#' Retrieves the interactions from one ligand-receptor resource
#' @inheritDotParams OmnipathR::import_post_translational_interactions
#' @inheritParams get_partners
#'
#' @importFrom OmnipathR import_post_translational_interactions
intercell_connections <- function(resource, ...){

    if(resource == 'OmniPath'){

        return(omnipath_intercell(...))

    }

    import_post_translational_interactions(
        resource = resource,
        entity_type = 'protein',
        ...
    ) %>%
        as_tibble() %>%
        mutate(category_intercell_source = "ligand",
               category_intercell_target = "receptor")

}


#' Retrieves ligands from one ligand receptor resource
#' @inheritDotParams intercell_connections
#' @inheritParams get_partners
get_ligands <- function(resource, ...){

    get_partners(side = 'ligand', resource = resource, ...)

}


#' Retrieves receptors from one ligand-receptor resource
#' @inheritDotParams intercell_connections
get_receptors <- function(resource, ...){

    get_partners(side = 'receptor', resource = resource, ...)

}



#' Retrieves intercellular communication partners (ligands or receptors) from
#' one ligand-receptor resource.
#' @inheritParams omnipath_partners
#' @param resource Name of current resource (taken from get_lr_resources)
#' @inheritDotParams omnipath_intercell
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
#' @param ... Passed to \code{OmnipathR::import_omnipath_intercell}.
#'
#' @importFrom OmnipathR import_omnipath_intercell
#' @importFrom magrittr %>%
omnipath_partners <- function(side, ...){

    causality <- list(ligand = 'trans', receptor = 'rec')

    import_omnipath_intercell(
        causality = causality[[side]],
        scope = 'generic',
        source = 'composite',
        entity_type = 'protein',
        ...
    )

}


#' Ligand-receptor data for the descriptive part
#'
#' Ligands, receptors and connections from each resource in a nested list
#' of tibbles.
#'
#' @seealso \code{\link{compile_ligrec}}
compile_ligrec_descr <- function(){
    compile_ligrec(lr_pipeline = FALSE)
}


#' Shuffle OmniPath Intercell DB
#' @param op_resource Intrcell DB to shuffle
#' @param .seed Value for set.seed
#' @return A shuffled omnipath-formatted resource
#' @import BiRewire tibble
#' @export
shuffle_omnipath <- function(op_resource,
                             .seed = 1004){

    ### These packages could go to "Suggests" in DESCRIPTION
    ### because not all users want to install all the tools
    ### to run one of them. Functions from these packages
    ### should be referred by :: to avoid warnings
    # library(BiRewire)
    set.seed(.seed)

    # make a vector proportional to the number of consensus directions
    stimul_num <- round(mean(op_resource$consensus_stimulation)/mean(op_resource$consensus_direction) * 100000)
    directed_vector <- append(rep(1, stimul_num), rep(-1,100000 - stimul_num))

    op_prep <- op_resource %>%
        # filter(entity_type_intercell_source != "complex",
        #        entity_type_intercell_target != "complex") %>%
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
    op_dsg = birewire.induced.bipartite(op_prep,
                                        delimitators = list(negative = '-1',
                                                            positive = '1'))
    op_sif = birewire.build.dsg(op_dsg,
                                delimitators = list(negative = '-1',
                                                    positive = '1'))

    # Rewire dsg
    random_dsg = birewire.rewire.dsg(dsg = op_dsg)
    random_sif = birewire.build.dsg(random_dsg,
                                    delimitators = list(negative = '-1',
                                                        positive = '1'))
    # Jacard dsg
    message(str_glue("Jaccard index between random and original resource: ",
                     {birewire.similarity.dsg(op_dsg, random_dsg)}))

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
