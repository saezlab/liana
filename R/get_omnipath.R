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
            'CellTalkDB',
            'OmniPath'
        )
    )

}



#' Function to get intercell resources as in OmniPath
#'
#' @return A list with intercellular resources from OmniPath
#'
#' @importFrom OmnipathR import_intercell_network
#' @importFrom magrittr %>%
#' @importFrom methods new
#' @export
get_omni_resources <- function(){

    omni_list <- get_lr_resources()

    # Need to come back and expand on categories
    # Replace ligand and receptor with more descriptive ones when available?
    # Is cell-cell contact taken into account?
    # adhesion, ECM, etc
    setClass(
        "OmniCriteria",
        slots = list(
            interactions_param="list",
            transmitter_param="list",
            receiver_param="list"
        )
    )


    # Get a list with dataframes of omnipath resources
    omni_resources <- omni_list %>%
        map(function(x){
            message(x)
            if(x!="OmniPath"){
                x_obj = methods::new(
                    "OmniCriteria",
                    interactions_param = list(resource = x),
                    transmitter_param = list(
                        resource = x,
                        category = "ligand"
                    ),
                    receiver_param = list(
                        resource = x,
                        category = "receptor"
                    )
                )

                import_intercell_network(
                    interactions_param = x_obj@interactions_param,
                    transmitter_param = x_obj@transmitter_param,
                    receiver_param = x_obj@receiver_param)
            } else{
                import_intercell_network(
                    transmitter_param = list(category = "ligand"),
                    receiver_param = list(category = "receptor")
                )
            }
        }) %>%
        setNames(omni_list)

    return(omni_resources)

}


#' Function to get unfiltered intercell resources
#'
#' For each resource and OmniPath variant compiles tables of ligands,
#' receptors and interactions.
#'
#' @details calls on \code{omnipath_intercell}, \code{intercell_connections},
#'     \code{get_partners}, and \code{intercell_connections}
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map
#'
#' @param omni_variants bool whether to get different OmniPath variants (e.g.
#'     _full, based on ligrec resource quality quartile, or if only lig_rec)
#' @param lr_pipeline bool whether to format for lr_pipeline and remove
#'     duplicate LRs (mainly from composite OmniDB due to category
#'     (adhesion vs lr))
#'
#' @return A list of OmniPath resources formatted according to the method pipes
#' @export
compile_ligrec <- function(omni_variants = FALSE, lr_pipeline = TRUE){

    # A list of OmniPath variants to be returned
    if(omni_variants){
        omnipath_variants <- list(
            OmniPath_q50 = list(quality = .5),
            OmniPath_ligrec = list(ligrec = TRUE),
            OmniPath_ligrec_q50 = list(ligrec = TRUE, quality = .5)
        )
    } else{
        omnipath_variants <- list()
    }

    log_success('Compiling ligand-receptor database data.')

    omni_resources <-
        get_lr_resources() %>%
        map(
            function(resource){
                list(
                    ligands = get_ligands(resource),
                    receptors = get_receptors(resource),
                    connections = intercell_connections(resource)
                )
            }
        ) %>%
        setNames(get_lr_resources()) %>%
        c(
            map(
                omnipath_variants,
                function(args){
                    args %<>% c(list(resource = 'OmniPath'))
                    list(
                        ligands = do.call(get_ligands, args),
                        receptors = do.call(get_receptors, args),
                        connections = do.call(intercell_connections, args)
                    )
                }
            )
        ) %>%
        {
            if(lr_pipeline) reform_omni(.)
            else .
        }

    log_success('Finished compiling ligand-receptor database data.')

    return(omni_resources)
}


#' Ligand-receptor data for the descriptive part
#'
#' Ligands, receptors and connections from each resource in a nested list
#' of tibbles.
#'
#' @seealso \code{\link{compile_ligrec}}
compile_ligrec_descr <- function(){

    compile_ligrec(omni_variants = TRUE, lr_pipeline = FALSE)

}


#' Helper Function to Reformat Omni_resources for LR Pipeline
#'
#' @param omni_resources List: output from \code{compile_ligrec}.
#'
#' @importFrom purrr map pluck
#' @importFrom dplyr distinct_at
#' @importFrom magrittr %>%
#'
#' @return A list of OmniPath resources, including OmniPath composite DB,
#'     a reshuffled OmniPath, and a Default with NULL ( tool pipelines run
#'     using their default resource)
reform_omni <- function(omni_resources){

    map(
        omni_resources,
        function(x){
            x %>%
            pluck("connections") %>%
            distinct_at(
                .vars = c(
                    "source_genesymbol", # remove duplicate LRs
                    "target_genesymbol"
                ),
                .keep_all = TRUE
            )
        }
    ) %>%
    append(
        list(
            "Random" = shuffle_omnipath(.$connectomeDB2020),
            "Default" = NULL
        ),
        .
    )

}



#' Retrieves intercellular interactions from OmniPath
#'
#' @param quality Numeric: a quality threshold for OmniPath between 0 and 1.
#' @param ligrec Use only ligand-receptor interactions from OmniPath.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter rename group_by mutate summarize_all inner_join
#' @importFrom dplyr ungroup first
#' @importFrom OmnipathR import_post_translational_interactions
#'
#' @return A tibble with Intercell interactions from OmnIPath
omnipath_intercell <- function(
    quality = NULL,
    ligrec = FALSE
){

    intracell <- c('intracellular_intercellular_related', 'intracellular')

    transmitters <-
        omnipath_partners(
            side = 'ligand',
            quality = quality,
            ligrec = ligrec
        ) %>%
        filter(!parent %in% intracell) %>%
        rename(category_source = source)

    receivers <-
        omnipath_partners(
            side = 'receptor',
            quality = quality,
            ligrec = ligrec
        ) %>%
        filter(!parent %in% intracell) %>%
        rename(category_source = source)

    interactions <- import_post_translational_interactions()

    interactions %>%
        inner_join(
            transmitters,
            by = c('source' = 'uniprot')
        ) %>%
        group_by(
            category, parent, source, target
        ) %>%
        mutate(
            database = paste(database, collapse = ';')
        ) %>%
        summarize_all(first) %>%
        inner_join(
            receivers,
            by = c('target' = 'uniprot'),
            suffix = c('_intercell_source', '_intercell_target')
        ) %>%
        group_by(
            category_intercell_source,
            parent_intercell_source,
            source,
            target,
            category_intercell_target,
            parent_intercell_target
        ) %>%
        mutate(
            database_intercell_target = paste(
                database_intercell_target,
                collapse = ';'
            )
        ) %>%
        summarize_all(first) %>%
        ungroup()

}




#' Retrieves the interactions from one ligand-receptor resource
#'
#' @param resource Name of the resource. For a full list of resources
#'     see \code{get_lr_resources}.
#' @param ... Passed to \code{omnipath_intercell} or
#'     \code{OmnipathR::import_post_translational_interactions}.
#'
#' @importFrom OmnipathR import_post_translational_interactions
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
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
#'
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
#'
#' @param side Either "ligand" or "receptor".
#' @param resource Name of the resource. For a full list of resources
#'     see \code{get_lr_resources}.
#' @inheritParams omnipath_partners
#' @inheritDotParams omnipath_intercell
#'
#' @importFrom rlang sym syms !!
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct rename
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
#' @param side Either 'ligand' (trans), 'receptor' (rec) or 'both' (both
#'     short or long notation can be used).
#' @param quality Numeric: a quality threshold for OmniPath between 0 and 1.
#' @param ligrec Use only ligand-receptor interactions from OmniPath.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter ungroup
#' @importFrom tibble as_tibble
#' @importFrom OmnipathR import_omnipath_intercell
omnipath_partners <- function(
    side,
    quality = NULL,
    ligrec = FALSE
){

    causality <- list(ligand = 'trans', receptor = 'rec')

    import_omnipath_intercell(
        causality = causality[[side]],
        scope = 'generic',
        source = 'composite',
        entity_type = 'protein') %>%
        as_tibble() %>%
        {`if`(
            is.null(quality),
            .,
            group_by(., parent) %>%
                filter(
                    consensus_score >= quantile(consensus_score, quality)
                ) %>%
                ungroup()
        )} %>%
        {`if`(
            ligrec,
            filter(., parent == side),
            .
        )}
}