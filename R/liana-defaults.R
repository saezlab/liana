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
    squidpy.defaults = NULL){

    default_args <- list(
        "cellchat" = cellchat.defaults %<>%
            `%||%`(list(
                nboot = 100,
                exclude_anns = NULL,
                thresh = 1,
                assay = "RNA",
                .normalize = TRUE,
                .do_parallel = FALSE,
                .raw_use = TRUE
                )),

        'connectome' = connectome.defaults %<>%
        `%||%`(list(
            .spatial = FALSE,
            min.cells.per.ident = 1,
            p.values = TRUE,
            calculate.DOR = FALSE,
            assay = 'RNA',
            .format = TRUE
        )),

        'italk' = italk.defaults %<>%
                `%||%`(list(
                    assay = 'RNA',
                    .format = TRUE,
                    .DE = TRUE
                )),

        'natmi' = natmi.defaults %<>%
            `%||%`(list(
                omnidbs_path = "data/input/omnipath_NATMI",
                natmi_path = "NATMI/",
                em_path = "data/input/test_em.csv",
                ann_path = "data/input/test_metadata.csv",
                output_path = "data/output/NATMI_test",
                assay = "RNA",
                .format = TRUE,
                .write_data = TRUE,
                .seed = 1004,
                .num_cor = 4)),

        'sca' = sca.defaults %<>%
            `%||%`(list(
                assay = 'RNA',
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))),

        'squidpy' = squidpy.defaults %<>%
           `%||%`(list(
               python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
               cluster_key="seurat_annotations",
               n_perms=1000,
               threshold=0.01,
               seed=as.integer(1004)
           ))
    )
}
