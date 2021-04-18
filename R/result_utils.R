#' Get top hits
#' @param spec_list list of spec objects with ligrec results
#' @return A list of top hits per tool/tool_parameter
get_top_hits <- function(spec_list, n_ints=c(100, 500, 5000)){
    map(n_ints, function(.tn){
        names(spec_list) %>%
            map(function(method_name){

                parnams <- names(spec_list[[method_name]]@method_scores)

                map(parnams, function(parm){
                    spec_list[[method_name]]@method_results %>%
                        map(function(res){
                            res %>%
                                distinct() %>%
                                top_n(n=if_else(spec_list[[method_name]]@method_scores[[parm]],
                                                .tn,
                                                -.tn),
                                      wt=!!sym(parm)) %>%
                                as_tibble()
                        })
                }) %>%
                    {if(length(spec_list[[method_name]]@method_scores) > 1)
                        setNames(., str_glue("{method_name}_{parnams}"))
                        else setNames(., method_name)
                    }
            }) %>% setNames(names(spec_list)) %>%
            flatten() %>%
            setNames(str_replace_all(names(.),
                                     pattern = "_",
                                     replacement = "\\."))
    }) %>% setNames(str_glue("top_{n_ints}"))
}




#' Get binary activity frequencies per cell-cell pair
#' @param sig_list List of significant hits per Method-resource combination.
#'  Named list of methods with each element being a named list of resources.
#' @return A tibble of cell pair frequencies, based on binarized activity
get_binary_frequencies <- function(sig_list){
    sig_list %>%
        enframe() %>%
        unnest(value) %>%
        mutate(name = map(names(sig_list), # get combined method and resource names
                          function(m_name){
                              map(names(sig_list[[m_name]]),
                                  function(r_name){
                                      str_glue("{m_name}_{r_name}")
                                  })
                          }) %>% unlist()) %>%
        mutate(value = value %>%
                   map(function(res) res %>%
                           unite(source, target, col = "clust_pair") %>%
                           group_by(clust_pair) %>%
                           summarise(cn = n()) %>%
                           mutate(freq = cn / sum(cn)) %>%
                           select(clust_pair, freq)
                   )) %>%
        unnest(value)
}


#' Convert list with MethodSpecifics objects to a dataframe with ranked
#' cell_pair frequencies.
#' @param spec_list list with appropriately populated MethodSpecifics objects
#' @return a tibble with cell_pair frequencies represented by the ranked
#' normalized scores for each method by cell_pair
#' @details See format_rank_frequencies for details
get_rank_frequencies <- function(spec_list){
    names(spec_list) %>%
        map(function(method_name){

            parnams <- names(spec_list[[method_name]]@method_scores)

            map(parnams, function(parm){
                spec_list[[method_name]]@method_results %>%
                    map(function(res){
                        res %>%
                            format_rank_frequencies(score_col=sym(parm),
                                                    .desc_order = spec_list[[method_name]]@method_scores[[parm]])
                    })
            }) %>%
                {if(length(spec_list[[method_name]]@method_scores) > 1)
                    setNames(., str_glue("{method_name}_{parnams}"))
                    else setNames(., method_name)
                }
        }) %>% setNames(names(spec_list)) %>%
        flatten() %>%
        setNames(str_replace_all(names(.),
                                 pattern = "_",
                                 replacement = "\\.")) %>%
        reform_rank_frequencies()
}


#' Helper function to get cell pair activity from rank averages
#'
#' @param result Result for a specific Tool-resource combinations
#' @param score_col Score column name provided by the tool
#' @param .desc_order Whether the most significant hits are in descending order,
#' i.e. the highest is the most sig
#' @details Cell pair ranks are averaged, then converted to z-scores, which
#' are multiplied by -1 as we want the lowest average ranks to have the highest
#' z scores
format_rank_frequencies <- function(result, score_col, .desc_order = TRUE){
    result %>%
        as_tibble() %>%
        filter(source!=target) %>%
        select(source, target, ligand, receptor, !!score_col) %>%
        filter(!is.nan(!!rlang::sym(score_col))) %>%
        mutate(edge_rank := if_else(rep(.desc_order, nrow(.)),
                                    min_rank(desc(!!rlang::sym(score_col))),
                                    min_rank(!!rlang::sym(score_col))))  %>%
        unite(source, target, col = "clust_pair") %>%
        group_by(clust_pair) %>%
        summarise(avg_rank = (mean(edge_rank))) %>%
        mutate(freq = -1*scale(avg_rank)[, 1])
}



#' Helper function to convert list with all resources ranked to frequencies df
#'
#' @param frequencies_list list with all resources ranked to frequencies df
reform_rank_frequencies <- function(frequencies_list){

    # Combine all results into tool_resource list
    lnames <- map(names(frequencies_list), function(l_name){
        map(names(frequencies_list[[l_name]]), function(r_name){
            str_glue("{l_name}_{r_name}")
        })
    }) %>% unlist()

    freq_df <- frequencies_list %>%
        purrr::flatten() %>%
        setNames(lnames) %>%
        enframe() %>%
        unnest(value)

    return(freq_df)
}



#' S4 Class used to format benchmark output.
#' @name MethodSpecifics-class
#'
#' @field method_name
#' @field method_results
#' @field method_scores
#'
#' @exportClass MethodSpecifics
setClass("MethodSpecifics",
         slots=list(method_name="character",
                    method_results = "list",
                    method_scores="list"))



#' Helper function to get cell number per cell type
#'
#' @param seurat_path Path to Seurat object of interest
#'
#' @return A tibble with Cell_subtype and Cell Number columns
#' @import Seurat dplyr tibble
get_cellnum <- function(seurat_path){
    crc_form <- readRDS(seurat_path)
    crc_meta <- crc_form@meta.data
    rm(crc_form)
    gc()

    # Cell Numbers
    crc_meta %>%
        select(Cell_clusters, Cell_subtype) %>%
        group_by(Cell_subtype) %>%
        summarise(cell_occur = n()) %>%
        arrange(Cell_subtype)
}



# I) Active Cell types


#' Function to get Activity per Cell Type heatmap
#'
#' @param top_list A resource-tool list with top hits for each combination.
#'
#' @return Cell Type Activity Heatmap
#'
#' @import pheatmap tidyverse
#' @inheritDotParams pheatmap::pheatmap
get_activecell <- function(top_list, ...){
    # Split by source/target cell
    top_frac <- top_list %>%
        map(function(db){
            db %>%
                enframe(name = "resource", value = "results") %>%
                mutate(results = results %>% map(function(res) res %>%
                                                     select(source, target))) %>%
                unnest(results) %>%
                pivot_longer(cols = c(source, target),
                             names_to = "cat",
                             values_to = "cell") %>%
                group_by(resource) %>%
                mutate(total_count = n()) %>%
                ungroup() %>%
                group_by(resource, cat, cell) %>%
                mutate(cell_count = n()) %>%
                group_by(resource, cat, cell) %>%
                mutate(cell_fraq = cell_count/total_count) %>%
                distinct() %>%
                unite(cell, cat, col = "cell_cat") %>%
                pivot_wider(id_cols = resource,
                            names_from = cell_cat,
                            values_from = cell_fraq,
                            values_fill = 0)
        }) %>%
        enframe(name = "method", value = "results_resource") %>%
        unnest(results_resource) %>%
        unite(method, resource, col = "mr") %>%
        mutate_all(~ replace(., is.na(.), 0)) %>%
        select(-c("Stromal cells_source",
                  "Stromal cells_target",
                  "Epithelial cells_target")) # FIX NATMI


    # annotation groups (sequential vectors as in heatmap_binary_list)
    method_groups <- top_frac %>%
        separate(mr, into = c("method", "resource"), sep = "_") %>%
        pull(method)
    resource_groups <- top_frac %>%
        separate(mr, into = c("method", "resource"), sep = "_") %>%
        pull(resource)

    # data frame with column annotations.
    # with a column for resources and a column for methods
    annotations_df <- data.frame(Resource = resource_groups,
                                 Method = method_groups)  %>%
        mutate(rn = top_frac$mr) %>%
        column_to_rownames("rn")

    annotations_row <- data.frame(cell_cat = colnames(top_frac)[-1]) %>%
        separate(cell_cat, sep="_", into = c("Cell", "Category"), remove = FALSE) %>%
        column_to_rownames("cell_cat") %>%
        select(Category)

    # List with colors for each annotation.
    mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                     Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))),
                     Category = c("#E41A1C", "#377EB8"))
    names(mycolors$Resource) <- unique(resource_groups)
    names(mycolors$Method) <- unique(method_groups)
    names(mycolors$Category) <- unique(annotations_row$Category)

    lab_rows <- annotations_row %>%
        rownames_to_column("cellname") %>%
        separate(cellname, into = c("cell", "cat"), sep = "_") %>%
        pull(cell)

    cellfraq_heat <- pheatmap(top_frac %>%
                                  column_to_rownames("mr") %>%
                                  t(),
                              annotation_row = annotations_row,
                              annotation_col = annotations_df,
                              annotation_colors = mycolors,
                              display_numbers = FALSE,
                              silent = FALSE,
                              show_colnames = FALSE,
                              color = colorRampPalette(c("darkslategray2",
                                                         "violetred2"))(20),
                              labels_row = ,
                              fontsize = 15,
                              drop_levels = TRUE,
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              border_color = NA,
                              treeheight_row = 0,
                              treeheight_col = 100,
                              ...)
}
