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


# colors from brewer.pal
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
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap %T>%
    uniq_per_res %T>%
    ligand_receptor_upset(upset_args = upset_args) %>%
    ligrec_classes_all %>%
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
#' @importFrom rlang !!! exec
figure_path <- function(fname, ...){

    args <- list(...)
    outdir <- args$outdir
    args$outdir <- NULL
    fname %<>% {exec(sprintf, ., !!!args)}

    outdir %>%
    ensure_outdir %>%
    file.path(fname)

}

#' Decomplexify specific resources that contain complexes
#' @param ligrec Nested list with ligand-receptor resources statistics, as
#'    produced by \code{compile_ligrec}.
#'
#' @return `ligrec` contents but complexes were split into rows of subunits
#'    and each interaction between a complex and another protein is duplicated
#'    for each subunit
ligrec_decomplexify <- function(ligrec){
    complex_resources <- c("CellChatDB",
                           "CellPhoneDB",
                           "Baccin2019",
                           "ICELLNET"
    )

    cats <- list(transmitters = "uniprot",
                 receivers = "uniprot",
                 interactions = c("source",
                                  "target"))

    ligrec <-
        ligrec %>% map2(names(.),
                        function(res, resname) map2(cats, names(cats),
                                                    function(col, cat)
                                                        if(resname %in% complex_resources){
                                                            decomplexify(res[[cat]], col)
                                                        } else{ res[[cat]] }
                        )
        )
    return(ligrec)
}

#' Helper Function to 'decomplexify' a column or vector of columns in a resource
#' @param resource a ligrec resource
#' @param column column to separate and pivot long (e.g. genesymbol or uniprot)
#'
#' @return returns a longer tibble with complex subunits on seperate rows
decomplexify <- function(resource, column){
    column %>%
        map(function(col){
            sep_cols <- c(str_glue("col{rep(1:5)}"))
            resource <<- resource %>%
                separate(col,
                         into = sep_cols,
                         sep = "_",
                         extra = "drop",
                         fill = "right") %>%
                pivot_longer(cols = sep_cols,
                             values_to = col,
                             names_to = NULL) %>%
                tidyr::drop_na(col) %>%
                distinct()
        })
    return(resource)
}




#' Helper Function to Convert a List of Resources into a binarized DF
#' @param interaction_list list resources with interactions alone
#' @returns Returns a 1/0 dataframe where 1 is assigned to interactions
#'    present in a given resource, 0 to absent interactions
binarize_resources <- function(interaction_list){
    map(names(interaction_list), function(l_name){
        interaction_list[[l_name]] %>%
            select(source, target) %>%
            unite("interaction", source, target, sep="_") %>%
            mutate(!!l_name := 1)
    }) %>% reduce(., full_join, by = "interaction") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
        as.data.frame()
}


#' Jaccard and Shared Elements heatmaps
#' @param ligrec List of lists with ligand-receptor data, as produced by
#'     \code{\link{compile_ligrec_descr}}. (to be changed to ligrec following
#'     ligrec_overlap)
#' @importFrom logger log_success
ligrec_overheats <- function(ligrec){

    log_success('Plotting Pairwise Jaccard and Overlap.')

    # To be extended to transmitters and receivers
    ligrec_binary <- ligrec %>%
        map(function(res) pluck(res, "interactions") %>%
                distinct_at(.vars = c("target",
                                      "source"))
            ) %>%
        binarize_resources()

    jaccheat_save(jacc_pairwise(ligrec_binary),
                  figure_path("interactions_jaccard_heat.pdf"),
                  "Jaccard Index")


    overheat_save(interactions_shared(ligrec_binary),
                  figure_path("interactions_shared_heat.pdf"),
                  "% Present")

    return(ligrec)
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
#' @importFrom dplyr n_distinct recode bind_rows n
ligrec_overlap <- function(ligrec){

    log_success('Finding overlaps between resources.')

    grp_vars <- list(
        transmitters = syms('uniprot'),
        receivers = syms('uniprot'),
        interactions = syms(c('source', 'target'))
    )

    keys <- c('transmitters', 'receivers', 'interactions')

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
                unique = ifelse(omnipath, op_unique, n_resources) == 1
            ) %>%
            shorten_resources()
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

#' Function to Calculate Unique Interactions Per Resource by Category
#' @inheritParams summarize_overlaps
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate_all group_by summarise ungroup bind_rows
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom logger log_success
uniq_per_res <- function(ligrec_olap){

    log_success('Summarizing unique stats.')

    ligrec_olap %>%
        map(function(cat){
            cat %>%
                group_by(resource, unique) %>%
                summarise(unq_or_not = n()) %>%
                ungroup() %>%
                pivot_wider(id_cols = resource,
                            names_from = unique,
                            values_from = unq_or_not) %>%
                mutate_all(~ replace(., is.na(.), 0)) %>%
                bind_rows(summarise_all(filter(.), ~if(is.numeric(.)) sum(.) else "Total")) %>%
                bind_rows(summarise_all(filter(., !(resource %in% c("OmniPath", "Total"))),
                                        ~if(is.numeric(.)) sum(.) else "Total_excl_Omni")) %>%
                mutate(unq_perc = `TRUE`/( `FALSE` + `TRUE`) * 100) %>%
                bind_rows(summarise_all(filter(., !(resource %in% c("OmniPath", "Total", "Total_excl_Omni"))),
                                        ~if(is.numeric(.)) median(.) else "Median_excl_Omni")) %>%
                bind_rows(summarise_all(filter(., !(resource %in% c("Total"))),
                                        ~if(is.numeric(.)) median(.) else "Median"))
        }
        ) %>%
        bind_rows(.id = "category") %>%
        select(category, resource, unq_perc) %>%
        pivot_wider(names_from = category, id_cols = resource, values_from = unq_perc) %>%
        write.csv(figure_path("uniqes_per_resource.csv"))
}



#' Function to get the intersect of vectors within the same list
#' @param vector_list List of character vectors (i.e. interactions per resource)
#' @param .names names of the list elements
#' @importFrom purrr map
get_intersect <- function(vector_list, .names){
    seq(length(vector_list)) %>%
        map(function(i)
            map(seq(length(vector_list)), function(j){
                length(intersect(
                    vector_list[[i]],
                    vector_list[[j]]
                ))
            }) %>% setNames(.names)
        )
}


#' Total vs. unique barplot
#'
#' Creates a barplot with the unique, shared and total number of ligands,
#' receptors and interactions in each resource.
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
#' @importFrom logger log_success
total_unique_bar <- function(ligrec_olap){

    log_success('Drawing overlap barplots.')

    res_order <-
        ligrec_olap$interactions %>%
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


#' Upset plots of ligands, receptors and interactions
#'
#' @importFrom purrr cross2 map walk
#' @importFrom magrittr %>%
#' @importFrom rlang quos
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
                    args$label == 'interactions',
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
            scale.intersections = "identity",
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
#' @importFrom OmnipathR annotated_network import_omnipath_annotations
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

    }else if(resource == 'OP-L'){

        return(ligrec %>% localization_ligrec_classes)

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

    ligrec$interactions %<>%
        annotated_network(annot = annot, !!attr) %>%
        filter(!!attr_src == !!attr_tgt) %>%
        select(-!!attr_tgt) %>%
        rename(!!attr := !!attr_src)

    annot %<>% select(uniprot, !!attr)

    ligrec$transmitters %<>%
        left_join(annot, by = 'uniprot')

    ligrec$receivers %<>%
        left_join(annot, by = 'uniprot')

    ligrec %<>% map(largest_groups, !!attr, largest = largest)

    return(ligrec)

}


#' Ligand-receptor classification from HGNC
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#' @param largest Numeric: how many of the largest groups to use. If `NULL`
#'     all groups will be used.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr left_join filter select rename
#' @importFrom purrr map
hgnc_ligrec_classes <- function(ligrec, largest = 15){

    hgnc_lig <- hgnc_annot('ligand')
    hgnc_rec <- hgnc_annot('receptor')

    ligrec$interactions %<>%
        left_join(hgnc_lig, by = c('source' = 'uniprot')) %>%
        left_join(
            hgnc_rec,
            by = c('target' = 'uniprot'),
            suffix = c('_source', '_target')
        ) %>%
        filter(category_source == category_target) %>%
        select(-category_target) %>%
        rename(category = category_source)

    ligrec$transmitters %<>% left_join(hgnc_lig, by = 'uniprot')
    ligrec$receivers %<>% left_join(hgnc_rec, by = 'uniprot')

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


#' Ligand-receptor classification by localization
#'
#' Classifies ligands, receptors and their interactions by localization:
#' either plasmam mebrane transmembrane, plasma membrane peripheral or
#' secreted. A receptor typically is not secreted, but here we keep this
#' option just to see how many of them are eventually annotated as secreted.
#' Such annotations are not always wrong, as some receptors indeed have
#' secreted forms. The interactions classified by all possible combinations
#' of ligand and receptor localizations.
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#'
#' @importFrom OmnipathR import_omnipath_intercell
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter select distinct mutate recode distinct
#' @importFrom rlang exec !!!
localization_ligrec_classes <- function(ligrec){

    locations <- list(
        secreted = 'S',
        plasma_membrane_transmembrane = 'T',
        plasma_membrane_peripheral = 'P'
    )

    location_groups <- list(
        P_P = "Direct",
        P_T = "Direct",
        T_P = "Direct",
        T_T = "Direct",
        S_P = "Secreted",
        S_T = "Secreted",
        T_S = "Others",
        S_S = "Others",
        P_S = "Others"
    )

    annot <-
        import_omnipath_intercell(
            aspect = 'locational',
            parent = locations %>% names,
            consensus_percentile = 50
        ) %>%
        filter(source == 'composite') %>%
        select(uniprot, location = category) %>%
        mutate(
            location = exec(recode, .x = location, !!!locations)
        ) %>%
        distinct

    ligrec$interactions %<>%
        left_join(annot, by = c('source' = 'uniprot')) %>%
        left_join(annot, by = c('target' = 'uniprot'),
                  suffix = c('_source', '_target')
        ) %>%
        filter(
            !is.na(location_source) &
                !is.na(location_target)
        ) %>%
        unite(location_source, location_target,
              col="location", sep = "_") %>%
        mutate(
            location = exec(recode, .x = location, !!!location_groups)
        )

    ligrec$transmitters %<>%
        left_join(annot, by = 'uniprot')

    ligrec$receivers %<>%
        left_join(annot, by = 'uniprot')

    return(ligrec)

}


#' Ligand-receptor data stacked barplots with classifications
#'
#' Calls \code{\link{ligrec_classes}} with classification from
#' many annotation resources.
#'
#' @param ligrec List of tibbles with ligand-receptor data, as produced by
#'     \code{\link{ligrec_overlap}}.
#'
#' @return Returns the input data frame unchanged.
#'
#' @importFrom magrittr %T>% %>%
#' @importFrom stringr str_sub str_to_title
ligrec_classes_all <- function(ligrec){

    ligrec %T>%
    {log_success('Stacked barplots of classifications.')} %T>%
    ligrec_classes_bar_enrich('SignaLink_pathway', pathway, NULL) %T>%
    ligrec_classes_bar_enrich('SIGNOR', pathway, 15) %T>%
    ligrec_classes_bar_enrich('NetPath', pathway, 15) %T>%
    ligrec_classes_bar_enrich('CancerSEA', state, 15) %T>%
    ligrec_classes_bar_enrich(
        'MSigDB',
        geneset,
        15,
        filter_annot = collection == 'hallmark',
        label_annot = function(x){str_to_title(str_sub(x, 10))}
    ) %T>%
    ligrec_classes_bar_enrich('HGNC', category, 15) %T>%
    ligrec_classes_bar_enrich('OP-L', location) %T>%
    {log_success('Finished stacked barplots of classifications.')} %>%
    invisible

}


#' Ligand-receptor data stacked barplots with classification
#'
#' Assigns classes to ligands, receptors and their interactions and for each
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
ligrec_classes_bar_enrich <- function(
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
    walk2(
        names(.),
        classes_bar_perc,
        resource,
        !!attr
    ) %>%
    walk2(
        names(.),
        classes_enrich,
        resource,
        !!attr
    ) %>%
    invisible

}


#' Stacked barplot from a data frame of classified entities
#'
#' @param data A data frame with classified entities (ligands, receptors or
#'     interactions).
#' @param entity The name of the entity, to be included in the output file
#'     name and the y axis label.
#' @param resource The name of the resource, to be included in the output
#'     file name.
#' @param var Name of the classifying variable.
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
            values = .palette1,
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



#' Percentage Stacked barplot from a data frame of classified entities
#'
#' @param data A data frame with classified entities (ligands, receptors or
#'     interactions).
#' @param entity The name of the entity, to be included in the output file
#'     name and the y axis label.
#' @param resource The name of the resource, to be included in the output
#'     file name.
#' @param var Name of the classifying variable.
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
classes_bar_perc <- function(data, entity, resource, var){

    var <- enquo(var)

    legend_title <- sprintf(
        '%s (%s)',
        var %>% quo_text %>% str_to_title,
        resource %>% str_replace('_.*$', '')
    )

    path <- figure_path('classes_perc_%s_%s.pdf', entity, resource)

    data %<>%
        filter(!is.na(!!var)) %>%
        order_by_group_size(resource) %>%
        mutate(
            !!var := factor(
                !!var,
                levels = sort(unique(!!var)),
                ordered = TRUE
            )
        ) %>%
        dplyr::group_by(resource, !!var) %>%
        summarise(n = n()) %>%
        mutate(perc = n / sum(n)) %>%
        ungroup()

    mean_perc <- data %>%
        dplyr::group_by(!!var) %>%
        summarise(perc = mean(perc), .groups = "keep") %>%
        mutate(resource = "Average")

    data %<>% bind_rows(mean_perc)

    p <- ggplot(data, aes(x = resource, y = perc,
                          fill = factor(!!var))) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
            values = .palette1,
            guide = guide_legend(title = legend_title)
        )  +
        ylab(str_to_title(entity)) +
        xlab('Resources') +
        geom_text(
            aes(label = round(perc*100)),
            position = position_stack(vjust = 0.5), size = 2, color = "white") +
        scale_y_continuous(labels = scales::percent) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45,
                                       size = 7, hjust = 1),
            axis.text.y = element_text(vjust = 1,
                                       size = 6, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size=6),
            legend.text = element_text(size=5),
            strip.text.x = element_blank(),
            legend.key.size = unit(3, 'mm')
        )

    wide <- (data %>% pull(!!var) %>% as.character %>% nchar %>% max) > 30

    cairo_pdf(path, width = `if`(wide, 7, 5), height = 3, family = 'DINPro')

    print(p)

    dev.off()

}


#' Enrichment dot plot from a data frame of classified entities
#'
#' @param data A data frame with classified entities (ligands, receptors or
#'     interactions).
#' @param entity The name of the entity, to be included in the output file
#'     name and the y axis label.
#' @param resource The name of the resource, to be included in the output
#'     file name.
#' @param var Name of the classifying variable.
#' @param ... Passed to \code{\link{enrich2}}.
#'
#' @return Returns `NULL`.
#'
#' @importFrom rlang ensym !! quo_text
#' @importFrom magrittr %>%
#' @importFrom purrr discard
#' @importFrom stringr str_to_title
#' @importFrom viridis scale_fill_viridis
#' @importFrom ggplot2 ggplot aes geom_tile xlab ylab
#' @importFrom ggplot2 theme_bw element_text theme
#' @importFrom dplyr pull
classes_enrich <- function(data, entity, resource, var, ...){

    var <- ensym(var)

    path <- figure_path('enrich_heatmap_%s_%s.pdf', entity, resource)

    data %<>%
        enrich2(!!var, ...) %>%
        arrange(desc(abs(enrichment))) %>%
        filter(
            !!var %in% head(
                unique(
                    filter(., !is.infinite(enrichment)) %>% pull(!!var)
                ),
                n = 15
            )
        ) %>%
        cluster_for_heatmap(!!var, resource, enrichment) %>%
        replace_inf(enrichment) %>%
        shorten_resources()

    lim <- data %>% pull(enrichment) %>% abs %>% max

    p <- ggplot(data, aes(x = resource, y = !!var, fill = enrichment)) +
        geom_tile() +
        scale_fill_viridis(
            option = 'cividis',
            limits = c(-lim, lim),
            guide = guide_colorbar(
                title = sprintf('Enrichment\nof %s', entity)
            )
        ) +
        xlab('Resources') +
        ylab(
            sprintf(
                '%s (%s)',
                var %>% quo_text %>% str_to_title,
                sub('_.*', '', resource)
            )
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
        )

    height <-
        data %>% pull(!!var) %>% levels %>% length / 7 + 1.6
    width <-
        data %>% pull(!!var) %>% as.character %>% nchar %>% max %>%
        {`if`(. > 25, `if`(. > 40, 7.5, 6.5), 5.5)}

    cairo_pdf(path, width = width, height = height, family = 'DINPro')

    print(p)

    dev.off()

    invisible(NULL)

}


#' Pairwise enrichment between two factors
#'
#' From a data frame with two categorical variables performs Barnard or
#' Fisher tests between all combinations of the levels of the two variables.
#' In each test the contingency table looks like: in \code{var1} belongs to
#' level x or not vs. in \code{var2} belongs to level y or not. By default
#' Fisher test is used simply because Barnard tests take very very long to
#' compute.
#'
#' @param data A data frame with entities labeled with two categorical
#'     variables \code{var1} and \code{var2}.
#' @param var1 Name of the first categorical variable.
#' @param var2 Name of the second categorical variable. Because in
#'     ligand-receptor data it's always `resource` the default is `resource`.
#' @param p_adj_method Adjustment method for multiple testing p-value
#'     correction (see \code{stats::p.adjust}).
#' @param test_method Characted: either "barnard", "barnard2" or "fisher".
#'     "barnard" uses \code{Barnard::barnard.test} while "barnard2" uses
#'     \code{DescTools::BarnardTest}.
#' @param ... Passed to the function executing the Barnard test.
#'
#' @return A data frame with all possible combinations of the two categorical
#'     variables, p-values, adjusted p-values and odds ratios.
#'
#' @importFrom rlang ensym !! !!! ensym quo_text exec
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter group_by group_modify mutate last
#' @importFrom purrr cross_df map
#' @importFrom Barnard barnard.test
#' @importFrom DescTools BarnardTest
#' @importFrom RCurl merge.list
enrich2 <- function(
    data,
    var1,
    var2 = resource,
    p_adj_method = 'fdr',
    test_method = 'fisher',
    ...
){

    var1 <- ensym(var1)
    var2 <- ensym(var2)
    var1_str <- quo_text(var1)
    var2_str <- quo_text(var2)

    t_f <- c('TRUE', 'FALSE')

    data %<>% select(!!var1, !!var2)

    data %>%
    {map(names(.), function(x){unique(.[[x]])})} %>%
    setNames(names(data)) %>%
    cross_df() %>%
    filter(!is.na(!!var1) & !is.na(!!var2)) %>%
    group_by(!!var1, !!var2) %>%
    group_modify(
        function(.x, .y){
            f1 <- factor(data[[var1_str]] == .y[[var1_str]][1], levels = t_f)
            f2 <- factor(data[[var2_str]] == .y[[var2_str]][1], levels = t_f)
            result <- fisher.test(f1, f2)
            odds_ratio <- result$estimate
            if(test_method == 'barnard'){
                sink('NUL')
                result <- exec(barnard.test, !!!table(f1, f2), ...)
                sink()
            }else if(test_method == 'barnard2'){
                param <-
                    list(...) %>%
                    merge.list(list(fixed = NA, method = 'boschloo'))
                result <- exec(BarnardTest, f1, f2, !!!param)
            }
            tibble(pval = last(result$p.value), odds_ratio = odds_ratio)
        }
    ) %>%
    ungroup() %>%
    mutate(
        padj = p.adjust(pval, method = p_adj_method),
        enrichment = ifelse(
            odds_ratio < 1,
            -1 / odds_ratio,
            odds_ratio
        ) %>% unname
    )

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


#' In a data frame shortens the names of certain resources
#'
#' @importFrom dplyr mutate recode
#' @importFrom magrittr %>%
#' @importFrom rlang !!!
shorten_resources <- function(data){

    data %>%
    mutate(
        resource = exec(recode, .x = resource, !!!.resource_short)
    )

}


#' Replace infinite values with a value slightly lower or higher than the
#' extrema
#'
#' @importFrom magrittr %>%
#' @importFrom rlang ensym !! :=
#' @importFrom dplyr pull mutate
#' @importFrom purrr discard
replace_inf <- function(data, var, s = 1.1){

    var <- ensym(var)

    n_inf <-
        data %>% pull(!!var) %>% discard(is.infinite) %>% min %>% min(0) * s
    p_inf <-
        data %>% pull(!!var) %>% discard(is.infinite) %>% max %>% max(0) * s

    data %>%
    mutate(
        !!var := ifelse(
            !!var == -Inf,
            n_inf,
            ifelse(
                !!var == Inf,
                p_inf,
                !!var
            )
        )
    )

}


#' Rearranges a data frame before plotting a heatmap
#'
#' Applies hierarchical clustering on two variables and converts them to
#' factors with the order of their levels corresponding to the clusters.
#'
#' @param data A data frame used for the heatmap.
#' @param cat1 Name of the first categorical variable (one axis of the
#'     heatmap).
#' @param cat2 Name of the second categorical variable (another axis of the
#'     heatmap).
#' @param val Name of the numeric variable corresponding to the color scale
#'     of the heatmap.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang ensym !! :=
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr pull mutate
cluster_for_heatmap <- function(data, cat1, cat2, val){

    cat1 <- ensym(cat1)
    cat2 <- ensym(cat2)
    val <- ensym(val)

    m <-
        data %>%
        replace_inf(!!val) %>%
        pivot_wider(
            id_cols = !!cat1,
            names_from = !!cat2,
            values_from = !!val
        ) %>%
        {`rownames<-`(
            as.matrix(select(., -1)),
            pull(., 1)
        )}

    cat1_ord <-
        m %>% dist %>% hclust %>% `$`('order') %>% {`[`(rownames(m), .)}
    cat2_ord <-
        m %>% t %>% dist %>% hclust %>% `$`('order') %>% {`[`(colnames(m), .)}

    data %>%
    mutate(
        !!cat1 := factor(!!cat1, levels = cat1_ord, ordered = TRUE),
        !!cat2 := factor(!!cat2, levels = cat2_ord, ordered = TRUE)
    )

}


#' Origin of major location attributes of intercellular interactions
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang exec !!! set_names
#' @importFrom dplyr group_by mutate select filter ungroup bind_rows
#' @importFrom dplyr summarize_all inner_join slice first
#' @importFrom purrr map cross2
#' @importFrom readr write_tsv
#' @export
localization_examples <- function(){

    locations <- c(
        'secreted',
        'plasma_membrane_transmembrane',
        'plasma_membrane_peripheral'
    )

    suffix <- c('_transmitter', '_receiver')
    src_up <- c('source' = 'uniprot')
    tgt_up <- c('target' = 'uniprot')

    icn <-
        import_intercell_network(entity_types = 'protein') %>%
        group_by(source, target) %>%
        mutate(
            categories_transmitter = paste(
                unique(sprintf(
                    '%s:%s',
                    category_intercell_source,
                    database_intercell_source
                )),
                collapse = '|'
            ),
            categories_receiver = paste(
                unique(sprintf(
                    '%s:%s',
                    category_intercell_target,
                    database_intercell_target
                )),
                collapse = '|'
            )
        ) %>%
        summarize_all(first) %>%
        ungroup() %>%
        select(
            source,
            target,
            source_genesymbol,
            target_genesymbol,
            is_directed,
            is_stimulation,
            is_inhibition,
            sources,
            references,
            categories_transmitter,
            categories_receiver
        )

    loc <-
        locations %>%
        map(
            function(lo){
                args <- list(
                    aspect = 'locational',
                    scope = 'generic',
                    entity_types = 'protein'
                )
                args[[lo]] = TRUE
                import_omnipath_intercell %>%
                exec(!!!args) %>%
                filter(parent != 'intracellular') %>%
                select(uniprot, database, category)
            }
        ) %>%
        set_names(locations)

    up_annot <-
        bind_rows(
            uniprot_annot(keyword),
            uniprot_annot(location),
            uniprot_annot(topology)
        ) %>%
        group_by(uniprot) %>%
        mutate(annot = paste(annot, collapse = ',')) %>%
        summarize_all(first) %>%
        ungroup()

    locations %>%
    cross2(.x = ., .y = .) %>%
    map(
        function(lo){
            loc_transm <- loc[[lo[[1]]]]
            loc_rec <- loc[[lo[[2]]]]
            icn %>%
            inner_join(
                loc_transm,
                by = src_up
            ) %>%
            inner_join(
                loc_rec,
                by = tgt_up,
                suffix = suffix
            ) %>%
            group_by(source, target) %>%
            mutate(
                locations_transmitter = paste(
                    unique(sprintf(
                        '%s:%s',
                        category_transmitter,
                        database_transmitter
                    )),
                    collapse = '|'
                ),
                locations_receiver = paste(
                    unique(sprintf(
                        '%s:%s',
                        category_receiver,
                        database_receiver
                    )),
                    collapse = '|'
                )
            ) %>%
            select(
                -category_transmitter,
                -category_receiver,
                -database_transmitter,
                -database_receiver
            ) %>%
            summarize_all(first) %>%
            ungroup() %>%
            left_join(up_annot, by = src_up) %>%
            left_join(up_annot, by = tgt_up, suffix = suffix) %>%
            mutate(major_location = exec(sprintf, '%s-%s', !!!lo)) %>%
            slice(sample(1:n()))

        }
    ) %>%
    bind_rows %T>%
    write_tsv('localization_examples.tsv')

}


#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_text !!
#' @importFrom OmnipathR import_omnipath_annotations
#' @importFrom dplyr select group_by mutate summarize_all ungroup first
#' @noRd
uniprot_annot <- function(var){

    var <- enquo(var)
    var_str <- quo_text(var)

    import_omnipath_annotations(
        resources = sprintf('UniProt_%s', var_str),
        entity_types = 'protein',
        wide = TRUE
    ) %>%
    select(uniprot, !!var) %>%
    group_by(uniprot) %>%
    mutate(annot = paste(!!var, collapse = ',')) %>%
    select(uniprot, annot) %>%
    summarize_all(first) %>%
    ungroup()

}


#' Overlap between resources pairwise
#'
#' @param ligrec_binary Binarized data frame with interactions.
#'
#' @return Shared elements between Resource Y and X dataframe
#'
#' @importFrom dplyr mutate filter distinct group_by group_nest ungroup
#' @importFrom dplyr mutate_at select rowwise across
#' @importFrom tidyr pivot_longer unnest pivot_wider
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom rlang !!!
#' @importFrom purrr map
#' @importFrom tidyselect starts_with
#' @importFrom forcats fct_relevel
interactions_shared <- function(ligrec_binary){
    interacts_per_resource <- ligrec_binary %>%
        as_tibble() %>%
        mutate(across(!starts_with("interaction"),
                      ~ifelse(.==1,interaction,.))) %>%
        pivot_longer(-interaction,
                     names_to = "resource",
                     values_to = "interact") %>%
        filter(interact != 0) %>%
        distinct()

    intersects_per_resource <- interacts_per_resource %>%
        select(resource, interaction = interact) %>%
        group_by(resource) %>%
        group_nest() %>%
        mutate(interaction = data %>% map(function(i) i$interaction)) %>%
        mutate(intersect = interaction %>% get_intersect(resource)) %>%
        rowwise() %>%
        mutate(resource_len = length(interaction)) %>%
        ungroup() %>%
        unnest(intersect)


    shared_per_resource <- intersects_per_resource %>%
        mutate(resource2 = names(intersect)) %>%
        unnest(intersect) %>%
        select(resource, resource2, intersect, resource_len) %>%
        mutate(shared_prop = intersect/resource_len) %>%
        select(resource, resource2, shared_prop)


    mean_shared <- shared_per_resource %>%
        filter(resource != resource2) %>%
        group_by(resource2) %>%
        summarise(shared_prop = mean(shared_prop)) %>%
        mutate(resource = "Mean Shared")

    shared_per_resource <- shared_per_resource %>%
        bind_rows(mean_shared) %>%
        distinct() %>%
        pivot_wider(
            id_cols = resource,
            names_from = resource2,
            values_from = shared_prop
        ) %>%
        pivot_longer(-resource) %>%
        as.data.frame() %>%
        mutate_at(vars(resource, "name"),
                  list(~recode(., .x=!!!.resource_short))) %>%
        mutate_if(is.character, as.factor)  %>%
        mutate(resource = fct_relevel(resource, "Mean Shared", after = Inf)) %>%
        mutate(value = value * 100)

    return(shared_per_resource)
}


#' Pairwise Jaccard Index between resources
#'
#' @param ligrec_binary
#'
#' @return Long data frame with Jaccard indeces.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate_at recode select
#' @importFrom tibble rownames_to_column
jacc_pairwise <- function(ligrec_binary){

    jacc_mat <-
        ligrec_binary %>%
        select(-interaction) %>%
        t() %>%
        get_simil_dist(.,
                       sim_dist = "simil",
                       method = "Jaccard",
                       diag = TRUE) %>%
        as.matrix()
    diag(jacc_mat) <- 1

    jacc_df <- jacc_mat %>%
        as.data.frame() %>%
        rownames_to_column("resource")  %>%
        pivot_longer(-resource) %>%
        as.data.frame() %>%
        mutate_at(vars(resource, "name"),
                  list(~recode(., .x=!!!.resource_short)))

    return(jacc_df)

}


#' Save Jaccard Heatmaps
#'
#' @param df Data frame to be pivoted and used to plot the heatmap.
#' @param plotname Name of the plot to be saved.
#'
#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 theme xlab theme_minimal element_text element_blank
#' @importFrom ggplot2 guide_colorbar geom_text
#' @importFrom grDevices cairo_pdf dev.off
#' @importFrom viridis scale_fill_viridis
#' @details to be merged with overheat_save
jaccheat_save <- function(df, plotname, guide_title){
    p <- ggplot(data = df) +
        geom_tile(aes(name, resource, fill = value)) +
        scale_fill_viridis(
            option = 'cividis',
            guide = guide_colorbar(
                title = guide_title
            )
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(vjust = 1, angle = 45,
                                       size = 16, hjust = 1),
            axis.text.y = element_text(vjust = 1,
                                       size = 16, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size=16),
            strip.text.x = element_blank(),
            panel.spacing = unit(1.5, "lines")
        ) +
        xlab("Resource") +
        geom_text(aes(name, resource, label = round(value, digits = 3)),
                  color = "white", size = 5) +
        scale_y_discrete(limits=rev)

    cairo_pdf(plotname, width = 16, height = 9, family = 'DINPro')

    print(p)

    dev.off()
}



#' Save Overlap Heatmaps
#'
#' @param df Data frame to be pivoted and used to plot the heatmap.
#' @param plotname Name of the plot to be saved.
#'
#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 theme xlab theme_minimal element_text element_blank
#' @importFrom ggplot2 guide_colorbar geom_text
#' @importFrom grDevices cairo_pdf dev.off
#' @importFrom viridis scale_fill_viridis
overheat_save <- function(df, plotname, guide_title){
    p <- ggplot(data = df) +
        geom_tile(aes(name, resource, fill = value)) +
        scale_fill_viridis(
            option = 'cividis',
            guide = guide_colorbar(
                title = guide_title
            )
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(vjust = 1, angle = 325,
                                       size = 16, hjust = 1),
            axis.text.y = element_text(vjust = 1,
                                       size = 16, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size=14),
            legend.text = element_text(size=12),
            strip.text.x = element_blank(),
            panel.spacing = unit(1, "lines")
        ) +
        xlab("Resource") +
        shadowtext::geom_shadowtext(aes(name, resource,
                                        label = round(value, digits = 1)),
                                    color = "white", size = 5, bg.colour='gray25') +
        facet_grid(.~name, scales='free_x', space="free_x") +
        scale_x_discrete(position = "top") +
        scale_y_discrete(limits=rev)

    cairo_pdf(plotname, width = 16, height = 9, family = 'DINPro')

    print(p)

    dev.off()
}
