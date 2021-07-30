# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

require(SingleCellExperiment)


liana_res <- liana_wrap(seurat_object,
                        method = c('natmi', 'connectome', 'logfc',
                                   'squidpy', 'sca', 'cellchat'),
                        resource = "OmniPath",
                        cellchat.params=list(
                            .normalize=TRUE
                        ))

liana_agg <- liana_res %>%
    liana_aggregate()

liana_agg



lr_res <- liana_pipe(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath")[[1]])

get_connectome(lr_res=lr_res)
