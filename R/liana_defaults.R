#' Function to pass Default Arguments for each method
#'
#' @param liana_pipe.params list of Parameters passed to NATMI \code{\link{liana_pipe}}
#' @param liana_call.params list of Parameters passed to NATMI \code{\link{liana_call}}
#' @param cellchat.params list of Parameters passed to CellChat \code{\link{call_cellchat}}
#' @param squidpy.params list of Parameters passed to Squidpy \code{\link{call_squidpy}}
#' @param call_connectome.params list of Parameters passed to Connectome \code{\link{call_connectome}}
#' @param call_italk.params list of Parameters passed to iTALK \code{\link{call_italk}}
#' @param call_natmi.params list of Parameters passed to NATMI \code{\link{call_natmi}}
#' @param call_sca.params list of Parameters passed to SingleCellSignalR \code{\link{call_sca}}
#' @param assay Assay name passed to `call_italk`, `call_sca`, `call_cellchat`,
#'    and `call_connectome`
#' @param decomplexify specify whether complexes in the resource should be
#'   dissociated and taken into account
#' @param expr_prop minimum proportion of gene expression per cell type (0.2 by default).
#'  One should consider setting this to an appropriate value between 0 and 1,
#'  as an assumptions of these methods is that communication is coordinated at the cluster level.
#' @param seed random seed integer
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#'  end of x before the mean is computed. This is relevant only for the
#'  CellPhoneDB algorithm re-implementation in liana.
#'
#'
#' @details The default parameters for each method can also be overwritten by
#'  manually passing a list of parameters for the appropraite method
#'   \code{\link{liana_wrap}}
#'
#' Further, each `get_*` method will by default obtain the default params passed
#'    via \code{\link{liana_pipe}} and \code{\link{liana_call}}. This is done so that most steps
#'    required for the calculation of these methods are undertaken only once.
#'
#' @return A list of the default parameters for each method
#'
#' @export
liana_defaults <- function(
    assay = "RNA",
    assay.type = "logcounts",
    decomplexify = TRUE,
    expr_prop = 0.2,
    seed = 1004,
    trim = 0,
    parallelize = FALSE,
    workers = 8,
    cellchat.params = NULL,
    squidpy.params = NULL,
    call_sca.params = NULL,
    cellphonedb.params = NULL,
    permutation.params = NULL,
    liana_pipe.params = NULL,
    liana_call.params = NULL,
    call_natmi.params = NULL,
    call_connectome.params = NULL,
    call_italk.params = NULL){

    default_args <- list(
        "cellphonedb" = cellphonedb.params %<>%
        `%||%`(
            list(
                workers = workers,
                parallelize = parallelize
            )
        ),

        # liana_scores (passed to get_* functions)
        "permutation" = permutation.params %<>%
            `%||%`(
                list(
                    nperms = 1000,
                    parallelize = parallelize,
                    workers = workers,
                    seed=seed,
                    trim = trim
                    )
            ),


        "liana_pipe" = liana_pipe.params %<>%
            `%||%`(
                list(
                    decomplexify = decomplexify,
                    test.type = "wilcox",
                    pval.type = "all",
                    expr_prop = expr_prop,
                    trim = trim,
                    assay = assay,
                    assay.type = assay.type
                    )
                ),

        # this thing needs to be either completely remove or moved to liana_wrap
        "liana_call" = liana_call.params %<>%
            `%||%`(
                list(
                    complex_policy = "min0",
                    decomplexify = decomplexify # should always be true
                    )
            ),

        # call_* functions/pipes
        "cellchat" = cellchat.params %<>%
            `%||%`(list(
                nboot = 100,
                expr_prop = expr_prop,
                exclude_anns = NULL,
                thresh = 1,
                assay = assay,
                .normalize = FALSE,
                .do_parallel = FALSE,
                .raw_use = TRUE
                )),

        'squidpy' = squidpy.params %<>%
           `%||%`(list(
               cluster_key = NULL,
               n_perms = 1000,
               threshold = expr_prop,
               seed = as.integer(seed),
               assay.type = assay.type
           )),

        # external call_* functions
        'call_sca' = call_sca.params %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))
                ),

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
                expr_file = "test_em.csv",
                meta_file = "metadata.csv",
                output_dir = "NATMI_results",
                assay = assay,
                assay.type = "logcounts",
                num_cor = 4,
                .format = TRUE,
                .write_data = TRUE,
                .use_raw = FALSE,
                .seed = seed,
                .natmi_path = NULL,
                .delete_output = FALSE
                )),

        'call_italk' = call_italk.params %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                .DE = TRUE
            ))
    )
}
