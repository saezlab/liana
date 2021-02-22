#' Benchmark Robustness Wrapper
#' @param subsampling Vector of subsampling (e.g. 0.9, 0.7, 0.5)
#' @param lr_call Ligad Receptor Algorithm to call
#' @param seurat_object seurat object to subsample and use in the analyses
#' @inheritDotParams lr_call functions
bench_robust <- function(subsampling,
                         lr_call,
                         seurat_object,
                         ...){

    .args <- list(...)

    map(subsampling, function(ss){
        exec(.fn = lr_call,
             seurat_object = seurat_subsample(seurat_object, subsampling = ss),
             !!!.args)
    })
}


#' Function to subsample Seurat object to fraction of cells by type
#' @param seurat_object Processed Seurat object with sc data
#' @param subsampling Value used to subsample cell clusters
#' @param .idents name of the idents
#' @param .seed Integer used to set the RNG seed
seurat_subsample <- function(seurat_object,
                             subsampling = 0.8,
                             .idents = "seurat_clusters",
                             .seed=1004){

    set.seed(.seed)

    features_to_keep <- seurat_object@meta.data %>%
        rownames_to_column("barcode") %>%
        as_tibble() %>%
        group_by(.data[[.idents]]) %>%
        mutate(clust_num = round(n() * subsampling)) %>%
        group_split() %>%
        map(function(x){
            sample(x$barcode, x$clust_num)
        }) %>%
        unlist()

    seurat_sub <- seurat_object[,features_to_keep]
    message(str_glue("Subsampling: {subsampling}"))
    seurat_sub@project.name <- str_glue("subsample_{subsampling}")

    return(seurat_sub)
}


#' Helper function used to format robustness results
#' @param lr_res Unformatted robustness results
#' @param db_list list of the OmniPath resources/dbs used in the subsampling
#' @param subsampling Vector of values used in the subsampling
#' @param .res_order bool whether the results are ordered by resource and
#' subsampling (e.g. OmniPath 1, OmniPath 0.8, OmniPath 0.6)
#'
#' @return A formatted tidy tibble with DB names and subsampling
robust_format_res <- function(lr_res,
                              db_list,
                              subsampling,
                              .res_order = TRUE){

    if(.res_order){
        .resource <- str_glue("{rep(names(db_list), each = length(subsampling))}_subsamp_{rep(subsampling, times = length(db_list))}")
    } else{
        .resource <- str_glue("{rep(names(db_list), times = length(subsampling))}_subsamp_{rep(subsampling, each = length(db_list))}")
    }

    lr_res%>%
        purrr::flatten() %>%
        enframe(name = "resource", value="lr_res") %>%
        mutate(resource = .resource) %>%
        mutate(resource = str_replace(resource, "\\_", "xx")) %>%
        separate(resource, into = c("resource", "subsample"), sep = "xx")
}

#' Function used to calculate the ROC curve of each robustness run
#' @param lr_form_res formatted LR results returned by robust_format_res()
#' @param predictor_metric Name of the predictor metric used in the
#' corresponding LR algorithm
#' @param predictor_thresh Threshold for the predictor metric
#' (e.g. 0.05 for pvalue)
#' @param .rank bool: whether to rank statistic by percentage
#' @return A Formatted tibble with calculated Robustness ROC
robust_get_roc <- function(lr_form_res,
                           predictor_metric,
                           predictor_thresh,
                           .rank = FALSE){
    # Get "Ground truth"
    lr_ground <- lr_form_res %>%
        filter(subsample == "subsamp_1") %>%
        select(resource, lr_res) %>%
        deframe() %>%
        map(function(ground){
            ground %>%
                mutate(predictor = ifelse(rep(!.rank, nrow(ground)),
                                          .data[[predictor_metric]],
                                          percent_rank(.data[[predictor_metric]]))) %>%
                mutate(truth = if_else(predictor <= predictor_thresh &
                !is.na(predictor), 1, 0)) %>%
                unite(ligand, receptor, source, target, col = "interaction") %>%
                select(interaction, truth)
        }) %>%
        enframe(name = "resource", value = "truth")

    lr_ground %>%
        unnest(truth) %>%
        filter(truth == 1) %>%
        group_by(resource) %>%
        summarise(count_truth = n()) %>%
        print()

    # calculate lr_roc
    lr_roc <- left_join(lr_form_res, lr_ground, by = "resource") %>%
        rowwise() %>%
        mutate(roc = list(calc_curve(lr_res,
                                     truth,
                                     predictor_metric = predictor_metric)))

    return(lr_roc)
}

