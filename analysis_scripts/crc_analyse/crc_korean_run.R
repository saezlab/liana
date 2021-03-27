library(intercell)
sapply(list.files("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/R/", pattern = ".R", full.names = TRUE), source)
setwd("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/")

# Korean CRC data
# sparsify_to_seurat(counts_loc = "input/crc_data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",
#                    meta_loc = "input/crc_data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",
#                    save_loc = "input/crc_data/crc_korean.rds",
#                    .format = TRUE)

crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
    format_crc_meta()
saveRDS("input/crc_data/crc_korean_modified.rds")

korclxx <- crc_korean %>%
    pluck("meta.data") %>%
    group_by(Cell_subtype) %>%
    summarise(n())

