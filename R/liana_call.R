#' Function to obtain SingleCellSignalR-like scores
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`LRscore`) as calculated
#'    by SingleCellSignalR
get_sca <- function(lr_res,
                    ...){
    liana_call(
        lr_res = lr_res,
        method = "sca",
        ...
    )
}

#' Function to obtain scConnect-like interaction scores
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with interaction scores (`interaction_score`)
#'   as calculated by scConnect
get_scconnect <- function(lr_res,
                          ...){
    liana_call(
        lr_res = lr_res,
        method = "scconnect",
        ...
    )
}



#' Function to obtain connectome-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`weight_sc`) as calculated
#'    by Connectome
get_connectome <- function(lr_res,
                           ...){
    liana_call(
        lr_res = lr_res,
        method = "connectome",
        ...
        )
}



#' Function to obtain NATMI-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`edge_specificity`)
#'    as calculated by NATMI
get_natmi <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "natmi",
        ...
    )
}



#' Function to obtain logFC weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with a logFC metric (`logfc_comb`). `logfc_comb` is
#'    calculated as the product of the (1 vs the rest) log2FC for each ligand
#'    and receptor gene
get_logfc <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "logfc",
        ...
    )
}


#' Function to obtain CellPhoneDB-like scores
#'
#' @inheritParams liana_scores
#' @inheritParams cellphonedb_score
#' @inheritDotParams cellphonedb_score
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
get_cellphonedb <- function(lr_res,
                            ...){
    liana_call(
        lr_res = lr_res,
        method = "cellphonedb",
        ...
    )
}


#' Function to obtain CytoTalk-like scores
#'
#' @inheritParams liana_scores
#' @inheritParams cellphonedb_score
#' @inheritDotParams cellphonedb_score
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
get_cytotalk <- function(lr_res,
                         ...){

    liana_call(
        lr_res = lr_res,
        method = "cytotalk",
        ...
    )
}


#' Wrapper Function to obtain scores via liana_pipe
#'
#' @inheritParams liana_pipe
#' @inheritDotParams liana_pipe
#' @inheritParams liana_scores
#' @inheritParams recomplexify
#'
#' @export
#'
#' @return lr_res modified to be method-specific
liana_call <- function(method,
                       seurat_object,
                       op_resource,
                       lr_res,
                       decomplexify = TRUE,
                       complex_policy = 'min0',
                       ...){

    liana_scores(.score_specs()[[method]],
                 lr_res = lr_res,
                 complex_policy = complex_policy,
                 decomplexify = decomplexify,
                 ...)
}



#' Function to obtain different scoring schemes
#'
#' @param score_object score_object specific to the test obtained from score_specs
#' @param lr_res ligand-receptor DE results and other stats between clusters
#' @param decomplexify whether to dissociate complexes into subunits and hence
#'    and hence take complexes into account (decomplexify) or not
#' @inheritParams recomplexify
#' @param ... dot params passed to `*_score` functions
#'
#' @return lr_res modified to be method-specific
liana_scores <- function(score_object,
                         lr_res,
                         decomplexify,
                         complex_policy,
                         ...){
    lr_res %<>%
        select(ligand, receptor,
               ends_with("complex"),
               source, target,
               ends_with("prop"),
               !!score_object@columns)

    if(decomplexify){
        lr_res %<>%
            recomplexify(
                lr_res = .,
                columns = score_object@columns,
                complex_policy = complex_policy)
    }

    args <-
        append(
            list(lr_res = lr_res,
            score_col = score_object@method_score),
            list(...)
        )

    exec(score_object@score_fun, !!!args) %>%
        ungroup() %>%
        select(-ends_with("prop"))
}


