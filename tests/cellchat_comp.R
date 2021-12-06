# Generate New Output
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                              "input", "testdata.rds"))

# Fix Default
liana_all <- liana_wrap(seurat_object,
                        resource = liana::show_resources()[c(2:6)], # [-c(1:2)] all resources except Default
                        method = c('call_natmi', 'call_connectome', 'logfc',
                                   'cellchat', 'call_sca', 'squidpy'),
                        expr_prop=0,
                        cellchat.params=list(.normalize=TRUE))
saveRDS(liana_all, "data/output/liana_default_resources.RDS")
