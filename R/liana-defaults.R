#' Default Arguments Helper Function
liana_defaults <- function(
    cellchat.defaults = NULL,
    connectome.defaults = NULL,
    italk.defaults = NULL,
    natmi.defaults = NULL,
    sca.defaults = NULL,
    squidpy.defaults = NULL,
    ...){

    mod_args <- list(...)

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

    default_args %>%
        map2(names(.), function(method, method_name){
            method %>%
                map2(names(.), function(arg, arg_name){
                    if(!is.null(mod_args[[method_name]][arg_name])){
                        mod_args[[method_name]][[arg_name]]
                    } else{
                        arg
                    }
                })
        })
}

liana_defaults()[["squidpy"]][["python_path"]]


xd <- liana_defaults(list("squidpy" = list("python_path" = ":)")))


xd <- liana_defaults(list("squidpy" = list("python_path" = ":)")))


xd$cellchat$nboot


xd[["squidpy"]][["python_path"]] <- "xx"


mod_args <- list("squidpy" = list("python_path" = ":)"))


mod_args %>%
    map2(names(.), function(method, method_name){
        method %>%
            map2(names(.), function(arg, arg_name){
                default_args[[method_name]][[arg_name]] <- ":)"
            })
    })


# options('cellchat.defaults')[[1]] %>%
#     `%||%`(list(
#         nboot = 100,
#         exclude_anns = NULL,
#         thresh = 1,
#         assay = "RNA",
#         .normalize = TRUE,
#         .do_parallel = FALSE,
#         .raw_use = TRUE
#     )) %T>%
#     options(cellchat.defaults = .)
#
# options('connectome.defaults')[[1]] %>%
#     `%||%`(list(
#         .spatial = FALSE,
#         min.cells.per.ident = 1,
#         p.values = TRUE,
#         calculate.DOR = FALSE,
#         assay = 'RNA',
#         .format = TRUE
#     )) %T>%
#     options(connectome.defaults = .)
#
#
# options('italk.defaults')[[1]] %>%
#     `%||%`(list(
#         assay = 'RNA',
#         .format = TRUE,
#         .DE = TRUE
#     )) %T>%
#     options(italk.defaults = .)
#
# options('natmi.defaults')[[1]] %>%
#     `%||%`(list(
#         omnidbs_path = "data/input/omnipath_NATMI",
            # natmi_path = "NATMI/",
            # em_path = "data/input/test_em.csv",
            # ann_path = "data/input/test_metadata.csv",
            # output_path = "data/output/NATMI_test",
            # assay = "RNA",
            # .format = TRUE,
            # .write_data = TRUE,
            # .seed = 1004,
            # .num_cor = 4
#     )) %T>%
#     options(natmi.defaults = .)
#
#
# options('sca.defaults')[[1]] %>%
#     `%||%`(
#         list(
#             assay = 'RNA',
#             .format = TRUE,
#             s.score = 0,
#             logFC = log2(1.5)
#             )) %T>%
#     options(sca.defaults = .)
#
#
# options('squidpy.defaults')[[1]] %>%
#     `%||%`(list(
#         python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
#         cluster_key="seurat_annotations",
#         n_perms=1000,
#         threshold=0.01,
#         seed=as.integer(1004)
#     )) %T>%
#     options(squidpy.defaults = .)
