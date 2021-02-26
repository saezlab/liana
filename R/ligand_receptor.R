#!/usr/bin/env Rscript

#
#  This file is part of the `intercell` R package
#
#  Copyright
#  2021
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  Author(s): Daniel Dimitrov
#             Charlotte Boys
#             Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Git repo: https://github.com/saezlab/Cell_Cell_Investigation
#

library(dplyr)
library(magrittr)
library(tibble)
library(rlang)
library(purrr)
library(tidyr)
library(ggplot2)
library(OmnipahR)

omnipath_variants <- list(
    OmniPath_full = list(),
    OmniPath_q50 = list(quality = .5),
    OmniPath_ligrec = list(ligrec = TRUE),
    OmniPath_ligrec_q50 = list(ligrec = TRUE, quality = .5)
)

intercell_resources <- c(
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
    'CellChatDB',
    'CellTalkDB',
    'connectomeDB2020',
    'Wojtowicz2020',
    'talklr'
)


intercell_param <- list(
    transmitter_all = list(
        causality = 'transmitter'
    ),
    transmitter_q50 = list(
        causality = 'transmitter',
        quality = .5
    ),
    transmitter_ligrec_only = list(
        causality = 'transmitter',
        only = 'ligand'
    ),
    transmitter_ligrec_only_q50 = list(
        causality = 'transmitter',
        only = 'ligand',
        quality = .5
    ),
    transmitter_cellphonedb = list(
        causality = 'transmitter',
        resource = 'CellPhoneDB'
    ),
    transmitter_ramilowski = list(
        causality = 'transmitter',
        resource = 'Ramilowski2015'
    ),
    receiver_all = list(
        causality = 'receiver'
    ),
    receiver_q50 = list(
        causality = 'receiver',
        quality = .5
    ),
    receiver_ligrec_only = list(
        causality = 'receiver',
        only = 'receptor'
    ),
    receiver_ligrec_only_q50 = list(
        causality = 'receiver',
        only = 'receptor',
        quality = .5
    ),
    receiver_cellphonedb = list(
        causality = 'receiver',
        resource = 'CellPhoneDB'
    ),
    receiver_ramilowski = list(
        causality = 'receiver',
        resource = 'Ramilowski2015'
    )
)


#' Retrieves the interactions from one ligand-receptor resource
intercell_connections <- function(resource, ...){

    if(resource == 'OmniPath'){

        return(omnipath_intercell(...))

    }

    import_post_translational_interactions(
        resource = resource,
        entity_type = 'protein',
        ...
    ) %>%
    as_tibble

}


#' Retrieves ligands from one ligand receptor resource
ligands <- function(resource, ...){

    partners(side = 'ligand', resource = resource, ...)

}


#' Retrieves receptors from one ligand-receptor resource
receptors <- function(resource, ...){

    partners(side = 'receptor', resource = resource, ...)

}


#' Retrieves intercellular communication partners (ligands or receptors) from
#' one ligand-receptor resource.
partners <- function(side, resource, ...){

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
        entity_type = 'protein'
    ) %>%
    as_tibble() %>%
    {`if`(
        is.null(quality),
        .,
        group_by(., parent) %>%
        filter(
            consensus_score >=
            quantile(consensus_score, quality)
        ) %>%
        ungroup()
    )} %>%
    {`if`(
        ligrec,
        filter(., parent == side),
        .
    )}

}


#' Retrieves intercellular interactions from OmniPath
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


#' For each resource and OmniPath variant compiles tables of ligands,
#' receptors and interactions
compile_ligrec <- function(){

    intercell_resources %>%
    map(
        function(resource){
            list(
                ligands = ligands(resource),
                receptors = receptors(resource),
                connections = intercell_connections(resource)
            )
        }
    ) %>%
    setNames(intercell_resources) %>%
    c(map(
        omnipath_variants,
        function(args){
            args %<>% c(list(resource = 'OmniPath'))
            list(
                ligands = do.call(ligands, args),
                receptors = do.call(receptors, args),
                connections = do.call(intercell_connections, args)
            )
        }
    ))

}


ligrec_overlap <- function(ligrec){

    grp_vars <- list(
        ligands = syms('uniprot'),
        receptors = syms('uniprot'),
        connections = syms(c('source', 'target'))
    )

    keys <- c('ligands', 'receptors', 'connections')

    keys %>%
    map(
        function(key){
            ligrec %>%
            map2(
                names(.),
                function(data, resource){
                    data[[key]] %>%
                    select(!!!grp_vars[[key]]) %>%
                    distinct() %>%
                    mutate(resource = resource)
                }
            ) %>%
            do.call(bind_rows, .) %>%
            mutate(omnipath = startsWith(resource, 'OmniPath')) %>%
            group_by(!!!grp_vars[[key]], omnipath) %>%
            mutate(n_resources = n()) %>%
            ungroup() %>%
            group_by(!!!grp_vars[[key]]) %>%
            mutate(op_unique = n_distinct(omnipath)) %>%
            ungroup() %>%
            mutate(unique = ifelse(omnipath, op_unique, n_resources) == 1)
        }
    ) %>%
    setNames(keys)

}


summarize_overlaps <- function(ligrec_olap){

    ligrec_olap %>%
    map(
        function(data, key){
            data %>%
            group_by(resource) %>%
            mutate(total = n()) %>%
            ungroup() %>%
            group_by(resource, unique) %>%
            mutate(n_unique = ifelse(unique, n(), NA)) %>%
            summarize_all(first) %>%
            ungroup() %>%
            select(resource, omnipath, total, n_unique) %>%
            pivot_longer(c(total, n_unique)) %>%
            group_by(resource, name) %>%
            mutate(value = max(value, 0, na.rm = TRUE)) %>%
            summarize_all(first) %>%
            ungroup()
        }
    )

}


total_unique_bar <- function(ligrec_olap){

    res_order <- ligrec_olap$connections %>%
        filter(name == 'total') %>%
        arrange(value) %>%
        pull(resource)

    ligrec_olap %>%
    map2(
        names(.),
        function(data, label){

            data %<>%
                mutate(
                    resource = factor(
                        resource,
                        levels = res_order,
                        ordered = TRUE
                    ),
                    name = factor(
                        name,
                        levels = c('total', 'n_unique'),
                        ordered = TRUE
                    )
                )

            p <- ggplot(data, aes(y = resource, x = value, fill = name)) +
                geom_col() +
                ylab('Records') +
                xlab('Resources') +
                ggtitle(str_to_title(label)) +
                theme_bw()

            cairo_pdf(
                sprintf('size_overlap_%s.pdf', label),
                width = 5,
                height = 8,
                family = 'DINPro'
            )

            print(p)

            dev.off()

        }
    ) %>%
    invisible()

}