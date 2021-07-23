# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

require(SingleCellExperiment)


liana_res <- liana_wrap(seurat_object,
                        method = c('natmi', 'connectome', 'logfc',
                                   'squidpy', 'sca'),
                        resource = "OmniPath")

liana_agg <- liana_res %>%
    liana_aggregate()
