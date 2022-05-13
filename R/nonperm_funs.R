#' Function Used to Calculate the Connectome-like `weight_sc` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#'
#' @noRd
#'
#' @return lr_res with an added `weight_sc` column
connectome_score <- function(lr_res,
                             score_col,
                             ...){
    expr_prop <- list(...)[["expr_prop"]]

    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := mean(c(ligand.scaled, receptor.scaled))) %>%
        filter(receptor.prop >= expr_prop & ligand.prop >= expr_prop)
}


#' Function Used to Calculate the NATMI-like `edge_specificity` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @return lr_res with an added `edge_specificity` column
#'
#' @noRd
#'
#' @details In the original NATMI implementation NAs are filtered out, but
#' here replace NAs with 0s, as they are needed to account for complexes
natmi_score <- function(lr_res,
                        score_col,
                        ...){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := ((ligand.expr*(ligand.sum^-1))) *
                    ((receptor.expr*(receptor.sum^-1)))) %>%
        mutate( {{ score_col }} := tidyr::replace_na(.data[[score_col]], 0))
}


#' Function Used to Calculate the logFC products (by default)
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @noRd
#'
#' @return lr_res with an added `logfc_comb` column
logfc_score <- function(lr_res,
                        score_col,
                        ...){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := mean(c(ligand.log2FC, receptor.log2FC)))
}



#' Function Used to Calculate the SigneCellSignalR `LRscore` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @noRd
#'
#' @return lr_res with an added `LRscore` column
sca_score <- function(lr_res,
                      score_col,
                      ...){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} :=
                    (ligand.expr^(1/2) * receptor.expr^(1/2))/
                    (global_mean + ligand.expr^(1/2) * receptor.expr^(1/2))
        )
}


### To be finished ----
#' Correlation Coefficient For Interactions
#'
#' @param sce SingleCellExperiment Object
#' @param lr_res a tabble with LR results obtained in the process of liana_pipe
#'
#' @noRd
#'
#' @return
corr_score <- function(lr_res,
                       sce,
                       ...){

    # should filter to cell type A and B first? - this way it's not specific
    # should also remove genes that are not in lr_res to save time
    corr_pairs <- scran::correlatePairs(sce) %>%
        as_tibble()
    corr_pairs <- corr_pairs %>%
        select(gene1=gene2,
               gene2=gene1,
               everything()) %>%
        bind_rows(corr_pairs) %>%
        select(gene1,
               gene2,
               rho,
               corr.FDR=FDR)

    corr_score <- lr_res %>%
        left_join(
            corr_pairs,
            by=c("ligand"="gene1",
                 "receptor"="gene2")
        ) %>%
        distinct()
}
