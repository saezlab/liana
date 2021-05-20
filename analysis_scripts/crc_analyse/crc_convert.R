require(Seurat)
require(tidyverse)

setwd("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/")

# Korean CRC data
sparsify_to_seurat(counts_loc = "input/crc_data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",
                   meta_loc = "input/crc_data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",
                   save_loc = "input/crc_data/crc_korean.rds")


# # Belgian CRC data
# sparsify_to_seurat(counts_loc = "input/crc_data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt",
#                    meta_loc = "input/crc_data/GSE144735_processed_KUL3_CRC_10X_annotation.txt",
#                    save_loc = "input/crc_data/crc_belgian.rds")
#
#
#
# # Different Protocols Data
# sparsify_to_seurat(counts_loc = "input/crc_data/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt",
#                    meta_loc = "input/crc_data/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt",
#                    save_loc = "input/crc_data/crc_protocols.rds")
