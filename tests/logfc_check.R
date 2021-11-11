liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

res1 <- liana_wrap(seurat_object,
                   method = c('logfc'),
                   resource = c('OmniPath'))
res1

