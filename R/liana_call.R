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
#' @param method name of the method to be called
#' @inheritParams liana_pipe
#' @inheritDotParams liana_pipe
#' @inheritParams liana_scores
#' @inheritParams recomplexify
#'
#' @export
#'
#' @return lr_res modified to be method-specific
liana_call <- function(method,
                       lr_res,
                       ...){

    liana_scores(.score_specs()[[method]],
                 lr_res = lr_res,
                 ...)
}



#' Function to obtain different scoring schemes
#'
#' @param score_object score_object specific to the test obtained from score_specs
#' @param lr_res ligand-receptor DE results and other stats between clusters
#' @param ... dot params passed to `*_score` functions
#'
#' @return lr_res modified to be method-specific
liana_scores <- function(score_object,
                         lr_res,
                         ...){

    # join all columns
    all_columns <- c(score_object@columns, score_object@add_columns)
    # Remove NAs for methods that don't have additional columns
    all_columns <- as.character(na.omit(all_columns))


    lr_res %<>%
        select(ligand, receptor,
               ends_with("complex"),
               source, target,
               ends_with("prop"),
               !!all_columns)

    lr_res %<>%
        recomplexify(
            lr_res = .,
            columns = score_object@columns,
            add_columns = score_object@add_columns,
            ...)  %>%
        # Select only the relevant columns
        select(source, target,
               ligand.complex, ligand,
               receptor.complex, receptor,
               ends_with("prop"),
               !!all_columns)

    args <-
        append(
            list(lr_res = lr_res,
            score_col = score_object@method_score),
            list(...)
        )

    # Get expr prop from defaults/kwargs
    expr_prop <- list(...)[["expr_prop"]]

    exec(score_object@score_fun, !!!args) %>%
        ungroup() %>%
        select(everything(), ends_with("prop")) %>%
        filter(receptor.prop >= expr_prop & ligand.prop >= expr_prop) %>%
        # ensure that there are no duplicates (edge cases where multiple subunits
        # have the same expr. - note that we also include method score to ensure
        # that no information is being lost + there are no issues)
        distinct_at(c("source", "target",
                      "ligand.complex", "receptor.complex",
                      score_object@method_score), .keep_all = TRUE)
}


