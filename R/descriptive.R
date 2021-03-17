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


#' Ligand-receptor resources descriptive plots
#'
#' Generates a number of comparative and descriptive plots about
#' ligand-receptor resources.
#'
#' @importFrom rlang %||%
#' @importFrom magrittr %>% %<>%
descriptive_plots <- function(ligrec = NULL, outdir = NULL){

    options(intercell.fig_desc_dir = outdir)

    ligrec %<>% `%||%`(compile_ligrec())

    ligrec %>%
    ligrec_overlap %>%
    summarize_overlaps %>%
    total_unique_bar

}


#' Makes sure the output directory exists
#'
#' @importFrom rlang %||%
#' @importFrom magrittr %>% %T>%
ensure_outdir <- function(outdir = NULL){

    outdir %>%
    `%||%`(options('intercell.fig_desc_dir')[[1]]) %>%
    `%||%`(file.path('figures', 'descriptive')) %T>%
    options(intercell.fig_desc_dir = .) %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)

}


#' Creates a path for a figure output
#'
#' @importFrom magrittr %>% %<>%
figure_path <- function(fname, ...){

    args <- list(...)
    outdir <- args$outdir
    args$outdir <- NULL
    fname %<>% {do.call(sprintf, c(list(.), args))}

    outdir %>%
    ensure_outdir %>%
    file.path(fname)

}


#' Overlaps between ligand-receptor resources
#'
#' Identifies the unique and shared ligands, receptors and interactions
#'
#' @param ligrec Nested list with ligand-receptor resources statistics, as
#'    produced by \code{compile_ligrec}.
#'
#' @return The same nested list of data frames as the input (`ligrec`) with
#'     new columns added: unique, op_unique, omnipath, n_resources.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!!
#' @importFrom purrr map map2
#' @importFrom dplyr mutate select distinct group_by ungroup n_distinct
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


#' Summary of overlaps between ligand-receptor resources
#'
#' Counts the unique and shared ligands, receptors and interactions
#'
#' @param ligrec_olap Nested list with ligand-receptor resources statistics
#'    with the shared and unique entities identified. Produced by
#'    \code{ligrec_overlap}.
#'
#' @return Summarized counts of unique and shared elements for each resource.
#'
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate ungroup summarize_all select first
#' @importFrom tidyr pivot_longer
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


#' Total vs. unique barplot
#'
#' Creates a barplot with the unique, shared and total number of ligands,
#' receptors and connections in each resource.
#'
#' @param ligrec_olap Summarized ligand-receptor overlaps, as produced by
#'     \code{summarize_overlaps}.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map2
#' @importFrom dplyr mutate filter arrange pull
#' @importFrom ggplot2 ggplot aes geom_col ylab xlab ggtitle theme_bw
#' @importFrom stringr str_to_title
total_unique_bar <- function(ligrec_olap){

    res_order <-
        ligrec_olap$connections %>%
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
                figure_path('size_overlap_%s.pdf', label),
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