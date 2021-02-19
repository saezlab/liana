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


#' Helper function used to
#'
#' @param df
#' @param ground "Ground Truth" - i.e. results considered significant
#' when using the whole dataset (no subsampling)
#'
#' @return tidy data frame with meta information for each subsampling
prepare_for_roc = function(df, ground, predictor_metric) {
    res = df %>%
        unite(ligand, receptor, source, target, col = "interaction") %>%
        left_join(ground, .) %>%
        mutate(response = case_when(truth == 1 ~ 1,
                                    truth == 0 ~ 0)) %>%
        mutate(predictor = abs(.[[predictor_metric]])) %>%
        filter(!is.na(predictor)) %>%
        select(interaction, response, predictor) %>%
        mutate(response = as.factor(response))
}





#' Helper Function to subsample Seurat object to fraction of cells by type
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




#' This function takes the elements of the `activity` column and calculates
#'    precision-recall and ROC curves (depending on `curve`).
#' The `activity` column is populated with the output for each stat method and
#'    results from the `run_benchmark()` function. Each of the elements
#'    in `activity` are results from runs of the \link{decouple} wrapper.
#'
#' @param df
#' @param downsampling logical flag indicating if the number of Negatives
#'    should be downsampled to the number of Positives
#' @param times integer showing the number of downsampling
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#' @param seed An integer to set the RNG state for random number generation. Use
#'    NULL for random number generation.
#'
#' @return tidy data frame with precision, recall, auc, n, cp, cn and coverage
#'    in the case of PR curve; or sensitivity and specificity, auc, n, cp, cn
#'    and coverage in the case of ROC.
#' @import yardstick
calc_curve = function(df,
                      ground,
                      predictor_metric = "pvalue",
                      downsampling = FALSE,
                      times = 1000,
                      curve = "ROC",
                      seed = 1004){

    if(curve=="PR"){
        res_col_1 <- "precision"
        res_col_2 <- "recall"
        curve_fun = yardstick::pr_curve
        auc_fun = yardstick::pr_auc
    }
    else if(curve=="ROC"){
        res_col_1 <- "sensitivity"
        res_col_2 <- "specificity"
        curve_fun = yardstick::roc_curve
        auc_fun = yardstick::roc_auc
    }

    df = df %>%
        prepare_for_roc(., ground, predictor_metric)


    if (sum(which(df$response == 0)) == nrow(df)){
        return(as_tibble(NULL))
    }

    cn = df %>% filter(.data$response == 0)
    cp = df %>% filter(.data$response == 1)

    if (downsampling == TRUE) {
        num_tp = nrow(cp)

        res = map_df(seq(from=1, to=times, by=1), function(i) {
            set.seed(seed)
            df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
                bind_rows(cp)

            r_sub = df_sub %>%
                curve_fun(.data$response, .data$predictor)

            auc = df_sub %>%
                auc_fun(.data$response, .data$predictor) %>%
                pull(.data$.estimate)

            res_sub = tibble({{ res_col_1 }} := r_sub %>% pull(res_col_1),
                             {{ res_col_2 }} := r_sub %>% pull(res_col_2),
                             th = r_sub$.threshold,
                             auc = auc,
                             n = length(which(df$response == 1)),
                             cp = nrow(cp),
                             cn = nrow(cn),
                             coverage = feature_coverage) %>%
                mutate("run" = i)

        })
        # Get Average AUC
        res$auc <- sum(res$auc)/length(res$auc)
        res$cn <- nrow(cp)

    } else {
        r = df %>%
            curve_fun(.data$response, .data$predictor)
        auc = df %>%
            auc_fun(.data$response, .data$predictor)

        res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
                     {{ res_col_2 }} := r %>% pull(res_col_2),
                     th = r$.threshold,
                     auc = auc$.estimate,
                     n = length(which(df$response == 1)),
                     cp = nrow(cp),
                     cn = nrow(cn)) %>%
            arrange(!!res_col_1, !!res_col_2)
    }

    return(res)
}
