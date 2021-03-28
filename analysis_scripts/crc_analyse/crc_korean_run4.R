library(intercell)
sapply(list.files("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/R/", pattern = ".R", full.names = TRUE), source)
setwd("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/")


# Load Data
crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
    format_crc_meta()


# Get Full Omni Resources
# omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")


# 6. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = crc_korean,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = TRUE,
                                   .do_parallel = FALSE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/crc_res/cellchat_results.rds")
