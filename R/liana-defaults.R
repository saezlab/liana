#' Function to pass Default Arguments for each method
#'
#' @param cellchat.params CellChat Parameters \link{liana::call_cellchat}
#' @param connectome.params Connectome Parameters \link{liana::call_connectome}
#' @param italk.params iTALK Parameters \link{liana::call_italk}
#' @param natmi.params NATMI Parameters \link{liana::call_natmi}
#' @param sca.params SingleCellSignalR Parameters \link{liana::call_sca}
#' @param squidpy.params Squidpy Parameters \link{liana::call_squidpy}
#'
#' @details The default parameters for each method can also be overwritten by
#'  manually passing a list of parameters for the appropraite method
#'   \link{liana::liana_wrap}
#'
#' @return A list of the default parameters for each method
#'
#' @export
liana_defaults <- function(
    cellchat.params = NULL,
    connectome.params = NULL,
    italk.params = NULL,
    natmi.params = NULL,
    sca.params = NULL,
    squidpy.params = NULL,
    assay = "RNA"){

    default_args <- list(
        "cellchat" = cellchat.params %<>%
            `%||%`(list(
                nboot = 1000,
                exclude_anns = NULL,
                thresh = 1,
                assay = assay,
                .normalize = FALSE,
                .do_parallel = FALSE,
                .raw_use = TRUE
                )),

        'connectome' = connectome.params %<>%
        `%||%`(list(
            min.cells.per.ident = 1,
            p.values = TRUE,
            calculate.DOR = FALSE,
            assay = assay,
            .format = TRUE,
            .spatial = FALSE
        )),

        'italk' = italk.params %<>%
                `%||%`(list(
                    assay = assay,
                    .format = TRUE,
                    .DE = TRUE
                )),

        'natmi' = natmi.params %<>%
            `%||%`(list(
                expr_file = "em.csv",
                meta_file = "metadata.csv",
                output_dir = "NATMI_test",
                assay = "RNA",
                num_cor = 4,
                .format = TRUE,
                .write_data = TRUE,
                .seed = 1004,
                .natmi_path = NULL
                )),

        'sca' = sca.params %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))),

        'squidpy' = squidpy.params %<>%
           `%||%`(list(
               cluster_key="seurat_annotations",
               n_perms=1000,
               threshold=0.01,
               seed=as.integer(1004)
           ))
    )
}



