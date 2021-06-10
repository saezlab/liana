#' Function to pass Default Arguments for each method
#'
#' @param cellchat.defaults CellChat Default Parameters \link{liana::call_cellchat}
#' @param connectome.defaults Connectome Default Parameters \link{liana::call_connectome}
#' @param italk.defaults iTALK Default Parameters \link{liana::call_italk}
#' @param natmi.defaults NATMI Default Parameters \link{liana::call_natmi}
#' @param sca.defaults SingleCellSignalR Default Parameters \link{liana::call_sca}
#' @param squidpy.defaults Squidpy Default Parameters \link{liana::call_squidpy}
#'
#' @details Default parameters for each method can also be passed manually
#'
#' @return A list of the default parametres for each method
#'
#' @export
liana_defaults <- function(
    cellchat.defaults = NULL,
    connectome.defaults = NULL,
    italk.defaults = NULL,
    natmi.defaults = NULL,
    sca.defaults = NULL,
    squidpy.defaults = NULL,
    assay = "RNA"){

    default_args <- list(
        "cellchat" = cellchat.defaults %<>%
            `%||%`(list(
                nboot = 1000,
                exclude_anns = NULL,
                thresh = 1,
                assay = assay,
                .normalize = FALSE,
                .do_parallel = FALSE,
                .raw_use = TRUE
                )),

        'connectome' = connectome.defaults %<>%
        `%||%`(list(
            min.cells.per.ident = 10,
            p.values = TRUE,
            calculate.DOR = FALSE,
            assay = assay,
            .format = TRUE,
            .spatial = FALSE
        )),

        'italk' = italk.defaults %<>%
                `%||%`(list(
                    assay = assay,
                    .format = TRUE,
                    .DE = TRUE
                )),

        'natmi' = natmi.defaults %<>%
            `%||%`(list(
                omnidbs_dir = "omnipath_NATMI",
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

        'sca' = sca.defaults %<>%
            `%||%`(list(
                assay = assay,
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))),

        'squidpy' = squidpy.defaults %<>%
           `%||%`(list(
               cluster_key="seurat_annotations",
               n_perms=1000,
               threshold=0.01,
               seed=as.integer(1004)
           ))
    )
}



