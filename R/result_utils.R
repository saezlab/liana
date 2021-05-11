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
                            parm_order <-
                                spec_list[[method_name]]@method_scores[[parm]]

                            res %>%
                                distinct() %>%
                                top_enh(n=if_else(parm_order,
                                                .tn,
                                                -.tn),
                                      wt=parm) %>%
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


#' Helper Function to handle specific cases for CellChat and Squidpy
#' @inheritDotParams dplyr::top_n
#' @importFrom dplyr top_n
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @return Ordered tibble/df as from top_n
top_enh <- function(...){

    elipses <- list(...)
    elipses$wt <- sym(elipses$wt)

    # Filter according to p-values for Squidpy pvalue and CellChat prob
    if(elipses$wt == "prob"){
        elipses[[1]] <- elipses[[1]] %>%
            filter(pval <= 0.05)
    } else if(elipses$wt == "pval"){
        elipses[[1]] <- elipses[[1]] %>%
            filter(pval <= 0.05)
    } else if(elipses$wt == "pvalue"){
        elipses[[1]] <- elipses[[1]] %>%
            filter(pvalue <= 0.05)
    } else if (elipses$wt == "LRscore"){
        elipses[[1]] <- elipses[[1]] %>%
            filter(LRscore >= 0.5)
    } else if (elipses$wt == "weight_sc"){
        elipses[[1]] <- elipses[[1]] %>%
            filter(p_val_adj.rec <= 0.05) %>%
            filter(p_val_adj.lig <= 0.05)
    }
    return(do.call(top_n, elipses))
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
#' @param frequencies_list list with all resources ranked to frequencies df
#' @return rank frequencies df
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
#' @field method_name name of the method (e.g. CellChat)
#' @field method_results Named list of method-resource results
#' @field method_scores Named list of the measures provided by the method
#'  and whether they should be interpreted in descending order as the value
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


#' Helper Function to get Similarties and Distances from binary dfs
#' @importFrom proxy dist simil
#' @return a similarity/dissimilarity matrix
get_simil_dist <- function(sim_dist = "simil", ...){
    do.call(sim_dist,
            list(...))
}


#' Helper Function to get a binary top hits DF (for all method-resource combos)
#' @param sig_list list of significant hits per method-resource combo
get_binary_df <- function(sig_list){
    # get method and resource names combined
    lnames <- map(names(sig_list), function(m_name){
        map(names(sig_list[[m_name]]), function(r_name){
            str_glue("{m_name}_{r_name}")
        })
    }) %>%
        unlist()

    # get binarized significant hits list (1 for sig per method, 0 if absent)
    binary_df <- sig_list %>%
        purrr::flatten() %>%
        setNames(lnames) %>%
        prepForUpset() %>%
        as_tibble() %>%
        column_to_rownames("interaction")
}



#' Get (Dis)Similarities per Resource/Method Combinations
#' @param sig_list list of top hits
#' @inheritDotParams proxy::simil
#' @import tibble
#' @import purrr
simdist_resmet <- function(sig_list,
                           ...){
    binary_df <- get_binary_df(sig_list)
    excl_res <- c(
        "Reshuffled",
        "Default"
        )

    methods <- names(sig_list)
    resources <- names(sig_list[[1]])
    resources <- resources[!(resources %in% excl_res)] # remove these

    # Get Sim/Diss between Methods
    method_sim <- methods %>% map(function(met){
        binary_df %>%
            select(starts_with(met)) %>%
            select(!ends_with(excl_res))
    }) %>% map(function(met_binary)
        get_simil_dist(
            x = t(met_binary),
            ...)) %>%
        setNames(methods)

    # Get Sim/Diss between Resources
    resource_sim <- resources %>% map(function(resource){
        binary_df %>%
            select(ends_with(resource))
    }) %>% map(function(res_binary)
        get_simil_dist(
            x = t(res_binary),
            ...)) %>%
        setNames(resources)

    # This can be extended with other combinations of (dis)similarity lists

    return(list(
        "meth" = method_sim,
        "reso"= resource_sim))
}



#' Get Similarity/Dissimilarity Stats from lists with binary matrices
#' @param ... Any list with with DFs (typically binarized) for which we wish to
#'  to calculate the mean, median, sd, and length.
#' @return a summary tibble
#' @import tibble
list_stats <- function(...){
    args <- list(...)

    # combined (i.e. vectorised simdists)
    df_comb <- args %>%
        enframe(value="simdist") %>%
        mutate(name = str_glue("{name}_comb")) %>%
        mutate(simdist = simdist %>%
                   map(function(sim_mat) as.vector(sim_mat) %>% unlist))
    # averaged simdists per each element in the list
    df_mean <- args %>%
        enframe(value="simdist") %>%
        mutate(name = str_glue("{name}_mean")) %>%
        rowwise() %>%
        mutate(simdist = list(simdist %>%
                                  map(function(sim_mat)
                                      mean(sim_mat))
                              %>% unlist)
        ) %>%
        ungroup()

    # bind and calculate averages
    df_stats <- bind_rows(df_comb, df_mean) %>%
        rowwise() %>%
        mutate(mn = mean(simdist),
               med = median(simdist),
               len = length(simdist),
               sd = sd(simdist),
               .min = min(simdist),
               .max = max(simdist)) %>%
        ungroup()

    return(df_stats)
}
