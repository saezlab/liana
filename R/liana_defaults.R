#' Function to pass Default Arguments for each method
#'
#' @param natmi.params \code{\link{get_natmi}}
#' @param connectome.params \code{\link{get_connectome}}
#' @param logfc.params \code{\link{get_logfc}}
#' @param cellchat.params CellChat Parameters \code{\link{call_cellchat}}
#' @param sca.params SingleCellSignalR Parameters \code{\link{call_sca}}
#' @param squidpy.params Squidpy Parameters \code{\link{call_squidpy}}
#' @param call_connectome.params Connectome Parameters \code{\link{call_connectome}}
#' @param call_italk.params iTALK Parameters \code{\link{call_italk}}
#' @param call_natmi.params NATMI Parameters \code{\link{call_natmi}}
#'
#' @details The default parameters for each method can also be overwritten by
#'  manually passing a list of parameters for the appropraite method
#'   \code{\link{liana_wrap}}
#'
#' @return A list of the default parameters for each method
#'
#' @export
liana_defaults <- function(
    assay = "RNA",
    cellchat.params = NULL,
    connectome.params = NULL,
    natmi.params = NULL,
    logfc.params = NULL,
    sca.params = NULL,
    squidpy.params = NULL,
    call_natmi.params = NULL,
    call_connectome.params = NULL,
    call_italk.params = NULL){

    default_args <- list(
        # liana_scores
        "connectome" = connectome.params %<>%
            `%||%`(list(

            )),

        "natmi" = natmi.params %<>%
            `%||%`(list(

            )),

        "logfc" = logfc.params %<>%
            `%||%`(list(

            )),


        # pipes
        "cellchat" = cellchat.params %<>%
            `%||%`(list(
                nboot = 100,
                exclude_anns = NULL,
                thresh = 1,
                assay = assay,
                .normalize = FALSE,
                .do_parallel = FALSE,
                .raw_use = TRUE
                )),

        'sca' = sca.params %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))),

        'squidpy' = squidpy.params %<>%
           `%||%`(list(
               cluster_key=NULL,
               n_perms=1000,
               threshold=0.01,
               seed=as.integer(1004)
           )),

        # deprecated call_* functions
        'call_connectome' = call_connectome.params %<>%
            `%||%`(list(
                min.cells.per.ident = 1,
                p.values = TRUE,
                calculate.DOR = FALSE,
                assay = assay,
                .format = TRUE,
                .spatial = FALSE
                )),

        'call_natmi' = call_natmi.params %<>%
            `%||%`(list(
                expr_file = "em.csv",
                meta_file = "metadata.csv",
                output_dir = "NATMI_results",
                assay = assay,
                num_cor = 4,
                .format = TRUE,
                .write_data = TRUE,
                .seed = 1004,
                .natmi_path = NULL
                )),

        'call_italk' = call_italk.params %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                .DE = TRUE
            ))
    )
}



