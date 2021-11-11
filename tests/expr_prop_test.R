# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
seurat_object <- readRDS("~/Repos/ligrec_decouple/data/input/citeseq/5k_pbmcs/5k_pbmcs_seurat.RDS")


# Run liana
liana_res <- liana_wrap(seurat_object,
                   method = c('call_connectome','connectome'),
                   resource = c('OmniPath'),
                   expr_prop = 0)
saveRDS(liana_res, "../ligrec_decouple/data/output/temp/connectome_test.RDS")

liana_conn <- liana_res$connectome %>%
    filter(weight_sc > 1.5)
og_conn <- liana_res$call_connectome %>%
    filter(weight_sc > 1.5)


# Expr prop filt
liana_conn_filt <- liana_wrap(seurat_object,
                              method = c('connectome'),
                              resource = c('OmniPath'),
                              expr_prop = 0.1)
liana_conn_filt_x <- liana_conn_filt %>%
    filter(weight_sc > 1.5)

saveRDS(liana_conn_filt, "../ligrec_decouple/data/output/temp/connectome_filt.RDS")








