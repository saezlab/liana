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
#' @importFrom magrittr %>% %<>% %T>%
descriptive_plots <- function(
    ligrec = NULL,
    outdir = NULL,
    upset_args = list()
){

    log_success('Descriptive visualization of ligand-receptor resources.')

    options(intercell.fig_desc_dir = outdir)

    ligrec %<>% `%||%`(compile_ligrec_descr())

    ligrec %>%
    ligrec_overlap %T>%
    ligand_receptor_upset(upset_args = upset_args) %>%
    summarize_overlaps %T>%
    total_unique_bar %T>%
    {log_success('Finished descriptive visualizations.')} %>%
    invisible


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
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang !!!
#' @importFrom purrr map map2
#' @importFrom dplyr mutate select distinct group_by ungroup n_distinct
ligrec_overlap <- function(ligrec){

    log_success('Finding overlaps between resources.')

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
    setNames(keys) %T>%
    {log_success('Finished finding overlaps between resources.')}

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

    log_success('Summarizing overlaps.')

    ligrec_olap %>%
    map(
        function(data){
            data %>%
            group_by(resource) %>%
            mutate(total = n()) %>%
            ungroup() %>%
            group_by(resource, unique) %>%
            mutate(
                n_unique = ifelse(unique, n(), NA),
                n_shared = ifelse(unique, NA, n())
            ) %>%
            summarize_all(first) %>%
            ungroup() %>%
            select(resource, omnipath, total, n_shared, n_unique) %>%
            pivot_longer(c(total, n_shared, n_unique)) %>%
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
#' @importFrom ggplot2 scale_fill_discrete guide_legend
#' @importFrom stringr str_to_title
total_unique_bar <- function(ligrec_olap){

    log_success('Drawing overlap barplots.')

    res_order <-
        ligrec_olap$connections %>%
        filter(name == 'total') %>%
        arrange(value) %>%
        pull(resource) %>%
        unique

    ligrec_olap %>%
    map2(
        names(.),
        function(data, label){

            data %<>%
                filter(name != 'total') %>%
                mutate(
                    resource = factor(
                        resource,
                        levels = res_order,
                        ordered = TRUE
                    ),
                    name = factor(
                        name,
                        levels = c('n_shared', 'n_unique'),
                        ordered = TRUE
                    )
                )

            p <- ggplot(data, aes(y = resource, x = value, fill = name)) +
                geom_col() +
                scale_fill_manual(
                    values = c('#B3C5E9', '#4268B3'),
                    label = c(n_shared = 'Shared', n_unique = 'Unique'),
                    guide = guide_legend(title = '')
                ) +
                xlab(str_to_title(label)) +
                ylab('Resources') +
                theme_bw()

            cairo_pdf(
                figure_path('size_overlap_%s.pdf', label),
                width = 5,
                height = 4,
                family = 'DINPro'
            )

            print(p)

            dev.off()

            return(data)

        }
    ) %>%
    invisible()

}


#' Upset plots of ligands, receptors and connections
#'
#' @importFrom purrr cross2 map walk
#' @importFrom magrittr %>%
ligand_receptor_upset <- function(data, upset_args = list()){

    log_success('Overlap upset plots.')

    data %>%
    names %>%
    cross2(c(TRUE, FALSE)) %>%
    map(setNames, c('label', 'omnipath')) %>%
    walk(
        function(args){
            args %>%
            c(
                list(
                    d = data[[args$label]],
                    upset_args = upset_args
                ),
                `if`(
                    args$label == 'connections',
                    quos(resource, source, target),
                    quos(resource, uniprot)
                )
            ) %>%
            do.call(what = upset_generic)
        }
    )

}


#' Creates an upset plot
#'
#' A wrapper to create an upset plot from a data frame
#'
#' @param data A data frame, one element from the output of
#'     \code{ligrec_overlap}.
#' @param label A label for the file name. Should refer to the entities we
#'     classify e.g. "ligand".
#' @param omnipath Logical: whether to include OmniPath.
#' @param upset_args List: additional arguments for \code{UpSetR::upset}.
#' @param ... Column names: the first one should be the set assignment, the
#'     rest define unique entities.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !! !!! enquos
#' @importFrom dplyr mutate select distinct filter group_by group_split
#' @importFrom dplyr group_keys pull
#' @importFrom purrr map
#' @importFrom UpSetR upset fromList
#' @importFrom RCurl merge.list
upset_generic <- function(data, label, omnipath, upset_args, ...){

    cols <- enquos(...)
    set_col <- cols[[1]]
    cols %<>% `[`(-1)

    data %<>%
    mutate(element = paste(!!!cols, sep = '__')) %>%
    select(!!set_col, element) %>%
    distinct() %>%
    {`if`(
        omnipath,
        filter(
            .,
            !startsWith(!!set_col, 'OmniPath') |
            !!set_col == 'OmniPath'
        ),
        filter(., !startsWith(!!set_col, 'OmniPath'))
    )} %>%
    group_by(!!set_col) %>%
    {setNames(
        group_split(., .keep = FALSE) %>% map(pull),
        group_keys(.) %>% pull
    )}

    upset_args %<>% merge.list(
        list(
            order.by = 'freq',
            mb.ratio = c(.5, .5),
            nsets = length(data),
            nintersects = 30,
            scale.intersections = `if`(omnipath, 'log10', 'identity'),
            show.numbers = 'no'
        )
    )

    path <- figure_path(
        'upset_%s_%s.pdf',
        label,
        `if`(omnipath, 'omnipath', 'no-omnipath')
    )

    cairo_pdf(path, width = 8, height = 4, family = 'DINPro')

        data %>%
        fromList %>%
        list %>%
        c(upset_args) %>%
        do.call(what = upset) %>%
        print()

    dev.off()

    invisible(data)

}


