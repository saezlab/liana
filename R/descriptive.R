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


# colors from Daniel
.palette1 <- c(
    '#8D4C6A',
    '#377EB8',
    '#419681',
    '#4DAF4A',
    '#727E76',
    '#984EA3',
    '#CB6651',
    '#FF7F00',
    '#FFBF19',
    '#FFFF33',
    '#D2AA2D',
    '#A65628',
    '#CE6B73',
    '#F781BF',
    '#E41A1C',
    '#1B9E77',
    '#D95F02',
    '#E7298A',
    '#66A61E',
    '#E6AB02'
)

# RWTH colors
.palette2 <- c(
    '#176FC1',
    '#A6D81C',
    '#F89D0E',
    '#00AAB0',
    '#0A6167',
    '#5B205F',
    '#7264B9',
    '#9E1639',
    '#D22027',
    '#ED0772',
    '#FFF200',
    '#4CBD38'
)

# 13 shades of 4 colors (calm blue)
.palette3 <- c(
    '#007E9D',
    '#00ACCB',
    '#9ED3DD',
    '#5F3950',
    '#A190B1',
    '#DCBAC7',
    '#2B4035',
    '#4E9495',
    '#7CA2A3',
    '#949673',
    '#BCB590',
    '#D7CBB5',
    '#6A5340'
)

# Paul Tol palette, distinct
.palette4 <- c(
    '#88CCEE',
    '#332288',
    '#44AA99',
    '#117733',
    '#999933',
    '#DDCC77',
    '#CC6677',
    '#882255',
    '#AA4499'
)

# Paul Tol palette, distinct 21 shades of 7 colors
.palette5 <- c(
    '#771155',
    '#AA4488',
    '#CC99BB',
    '#114477',
    '#4477AA',
    '#77AADD',
    '#117777',
    '#44AAAA',
    '#77CCCC',
    '#117744',
    '#44AA77',
    '#88CCAA',
    '#777711',
    '#AAAA44',
    '#DDDD77',
    '#774411',
    '#AA7744',
    '#DDAA77',
    '#771122',
    '#AA4455',
    '#DD7788'
)

.resource_short <- list(
    Ramilowski2015 = 'Ramilowski',
    connectomeDB2020 = 'ConnDB2020',
    Guide2Pharma = 'GuidePharm',
    OmniPath_ligrec_q50 = 'OmniPath_LRQ',
    OmniPath_ligrec = 'OmniPath_LR',
    OmniPath_q50 = 'OmniPath_Q'
)


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
    ligand_receptor_classes_bar_all %>%
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
#' @importFrom rlang !!! exec
#' @importFrom purrr map map2
#' @importFrom dplyr mutate select distinct group_by ungroup
#' @importFrom dplyr n_distinct recode
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
            mutate(
                unique = ifelse(omnipath, op_unique, n_resources) == 1,
                resource = exec(recode, .x = resource, !!!.resource_short)
            )
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
#' @importFrom grDevices cairo_pdf
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
#' @importFrom grDevices cairo_pdf
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
            show.numbers = 'no',
            sets.x.label = str_to_title(label),
            mainbar.y.label = sprintf('Shared %s', label)
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


#' Combines ligand-receptor and classification data
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#' @param resource Character: name of the annotation resource.
#' @param attr Name of the classifying variable from the annotation resource
#'     (e.g. pathway).
#' @param largest Numeric: how many of the largest groups to use. If `NULL`
#'     all groups will be used.
#' @param filter_annot Expression for filtering the annotation data frame.
#' @param label_annot Function to process the labels of the annotation
#'     variable (e.g. to shorten long strings).
#'
#' @importFrom rlang enquo quo_text sym !! :=
#' @importFrom OmnipathR annotated_network
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr rename filter select left_join
#' @importFrom purrr map
ligand_receptor_classes <- function(
    ligrec,
    resource,
    attr,
    largest = NULL,
    filter_annot = TRUE,
    label_annot = identity
){

    if(resource == 'HGNC'){

        return(ligrec %>% hgnc_ligrec_classes(largest = largest))

    }

    attr <- enquo(attr)
    attr_str <- quo_text(attr)
    attr_src <- sprintf('%s_source', attr_str) %>% sym
    attr_tgt <- sprintf('%s_target', attr_str) %>% sym
    filter_annot <- enexprs(filter_annot)

    annot <-
        import_omnipath_annotations(resource = resource, wide = TRUE) %>%
        filter(!!!filter_annot) %>%
        mutate(!!attr := label_annot(!!attr))

    ligrec$connections %<>%
        annotated_network(annot = annot, !!attr) %>%
        filter(!!attr_src == !!attr_tgt) %>%
        select(-!!attr_tgt) %>%
        rename(!!attr := !!attr_src)

    annot %<>% select(uniprot, !!attr)

    ligrec$ligands %<>%
        left_join(annot, by = 'uniprot')

    ligrec$receptors %<>%
        left_join(annot, by = 'uniprot')

    ligrec %<>% map(largest_groups, !!attr, largest = largest)

    return(ligrec)

}


#' Ligand-receptor classification from HGNC
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr left_join filter select rename
#' @importFrom purrr map
hgnc_ligrec_classes <- function(ligrec, largest = 15){

    hgnc_lig <- hgnc_annot('ligand')
    hgnc_rec <- hgnc_annot('receptor')

    ligrec$connections %<>%
        left_join(hgnc_lig, by = c('source' = 'uniprot')) %>%
        left_join(
            hgnc_rec,
            by = c('target' = 'uniprot'),
            suffix = c('_source', '_target')
        ) %>%
        filter(category_source == category_target) %>%
        select(-category_target) %>%
        rename(category = category_source)

    ligrec$ligands %<>% left_join(hgnc_lig, by = 'uniprot')
    ligrec$receptors %<>% left_join(hgnc_rec, by = 'uniprot')

    ligrec %<>% map(largest_groups, category, largest = largest)

    return(ligrec)

}


#' Classes from HGNC
#'
#' @importFrom OmnipathR import_omnipath_intercell
#' @importFrom magrittr %>%
#' @importFrom dplyr select
hgnc_annot <- function(parent){

    import_omnipath_intercell(
        resource = 'HGNC',
        scope = 'specific',
        parent = parent,
        entity_type = 'protein'
    ) %>%
    select(uniprot, category)

}


#' Ligand-receptor data stacked barplots with classifications
#'
#' Calls \code{\link{ligand_receptor_classes_bar}} with classification from
#' many annotation resources.
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#'
#' @return Returns the input data frame unchanged.
#'
#' @importFrom magrittr %T>% %>%
#' @importFrom stringr str_sub str_to_title
ligand_receptor_classes_bar_all <- function(ligrec){

    ligrec %T>%
    {log_success('Stacked barplots of classifications.')} %T>%
    ligand_receptor_classes_bar('SignaLink_pathway', pathway, NULL) %T>%
    ligand_receptor_classes_bar('SIGNOR', pathway, 15) %T>%
    ligand_receptor_classes_bar('NetPath', pathway, 15) %T>%
    ligand_receptor_classes_bar('CancerSEA', state, 15) %T>%
    ligand_receptor_classes_bar(
        'MSigDB',
        geneset,
        15,
        filter_annot = collection == 'hallmark',
        label_annot = function(x){str_to_title(str_sub(x, 10))}
    ) %T>%
    ligand_receptor_classes_bar('HGNC', category, 15) %T>%
    {log_success('Finished stacked barplots of classifications.')} %>%
    invisible

}


#' Ligand-receptor data stacked barplots with classification
#'
#' Assigns classes to ligands, receptors and their connections and for each
#' of these entities creates a stacked barplot.
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#' @param resource Character: name of the annotation resource.
#' @param attr Name of the classifying variable from the annotation resource
#'     (e.g. pathway).
#' @param largest Numeric: how many of the largest groups to use. If `NULL`
#'     all groups will be used.
#' @param ... Passed to \code{\link{ligand_receptor_classes}}.
#'
#' @importFrom rlang enquo !!
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr walk2
ligand_receptor_classes_bar <- function(
    ligrec,
    resource,
    attr,
    largest = NULL,
    ...
){

    attr <- enquo(attr)

    ligrec %<>%
    ligand_receptor_classes(
        resource,
        !!attr,
        largest = largest,
        ...
    ) %>%
    walk2(
        names(.),
        classes_bar,
        resource,
        !!attr
    ) %>%
    invisible

}


#' Stacked barplot from a data frame of classified entities
#'
#' @param data A data frame with classified entities (ligands, receptors or
#'     connections).
#' @param entity The name of the entity, to be included in the output file
#'     name and the y axis label.
#' @param resource The name of the resource, to be included in the output
#'     file name.
#'
#' @return Returns `NULL`.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang enquo !! quo_text :=
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab theme_bw theme
#' @importFrom ggplot2 scale_fill_manual guide_legend element_text
#' @importFrom dplyr filter mutate pull
#' @importFrom stringr str_to_title str_replace
#' @importFrom grDevices cairo_pdf
classes_bar <- function(data, entity, resource, var){

    var <- enquo(var)
    legend_title <- sprintf(
        '%s (%s)',
        var %>% quo_text %>% str_to_title,
        resource %>% str_replace('_.*$', '')
    )

    path <- figure_path('classes_%s_%s.pdf', entity, resource)

    data %<>%
        filter(!is.na(!!var)) %>%
        order_by_group_size(resource) %>%
        mutate(
            !!var := factor(
                !!var,
                levels = sort(unique(!!var)),
                ordered = TRUE
            )
        )

    p <- ggplot(data, aes(x = resource, fill = !!var)) +
        geom_bar() +
        stat_count() +
        scale_fill_manual(
            values = .palette5,
            guide = guide_legend(title = legend_title)
        ) +
        xlab('Resources') +
        ylab(str_to_title(entity)) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
            legend.text = element_text(size = 7),
            legend.key.size = unit(3, 'mm')
        )

    wide <- (data %>% pull(!!var) %>% as.character %>% nchar %>% max) > 30

    cairo_pdf(path, width = `if`(wide, 7, 5), height = 3, family = 'DINPro')

        print(p)

    dev.off()

}


#' Keep only the largest groups according to a grouping variable
#'
#' @param data A data frame.
#' @param var Name of the grouping variable.
#' @param largest Numeric: how many of the largest groups to keep. If `NULL`
#'     the data frame will be returned unchanged.
#'
#' @importFrom rlang ensym !!
#' @importFrom magrittr %>%
#' @importFrom dplyr add_count filter arrange desc filter select
#' @importFrom utils head
largest_groups <- function(data, var, largest = NULL){

    var <- ensym(var)

    data %>%
    {`if`(
        is.null(largest),
        .,
        add_count(., !!var) %T>%
        arrange(desc(n)) %>%
        filter(!!var %in% head(unique(!!var), n = largest)) %>%
        select(-n)
    )}

}


#' Orders the levels of a factor variable by their number of elements
#'
#' @param data A data frame.
#' @param var Name of a variable in the data frame.
#'
#' @return The data frame with the variable converted to a factor and its
#'     levels ordered.
#'
#' @importFrom rlang enquo !! :=
#' @importFrom dplyr add_count arrange desc mutate select
#' @importFrom magrittr %>%
order_by_group_size <- function(data, var){

    var <- enquo(var)

    data %>%
    add_count(resource) %>%
    arrange(desc(n)) %>%
    mutate(
        !!var := factor(
            !!var,
            levels = unique(!!var),
            ordered = TRUE
        )
    ) %>%
    select(-n)

}