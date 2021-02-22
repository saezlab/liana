#' Helper function used to format LR results to the appropriate format for the
#' calc_curve function
#'
#' @param lr_form_res Formatted LR results with subsampling and the resources
#' used for each tool.
#' @param ground "Ground Truth" - i.e. results considered significant
#' when using the whole dataset (no subsampling)
#'
#' @return A tidy data frame used to assess robustness
prepare_for_roc = function(lr_form_res, ground, predictor_metric) {
    res = lr_form_res %>%
        unite(ligand, receptor, source, target, col = "interaction") %>%
        left_join(ground, ., by = "interaction") %>%
        mutate(response = case_when(truth == 1 ~ 1,
                                    truth == 0 ~ 0)) %>%
        mutate(predictor = .[[predictor_metric]]) %>%
        filter(!is.na(predictor)) %>%
        select(interaction, response, predictor) %>%
        mutate(response = as.factor(response))
}




#' This function takes the elements of the lr_res column and calculates
#' robustness. The way that we assess robustness is by subsampling clusters.
#' Significant hits from the full data set are regarded as the ground truth
#' and results from subsampled data sets are evaluated against the ground truth.
#'
#' @param lr_form_res Formatted LR results with subsampling and the resources
#' used for each tool.
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


#' Get Robustness Roc Plot
#' @param lr_format_roc formatted LR results with ROC calculated
#' @return ggplot2 object ROC plot
robust_roc_plot <- function(lr_format_roc){
    lr_format_roc %>%
        unite(resource, subsample, col = "resource_subsample") %>%
        select(resource_subsample, roc) %>%
        unnest(roc) %>%
        ggplot(., aes(x = 1-specificity,
                               y = sensitivity,
                               colour = resource_subsample)) +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        xlab("FPR (1-specificity)") +
        ylab("TPR (sensitivity)")
}
