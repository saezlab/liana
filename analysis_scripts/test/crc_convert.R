#
crc_protocols_counts <- read_delim("input/crc_data/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt",
                         delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")


crc_korean_counts <- read_delim("input/crc_data/GSE132465_GEO_processed_head.txt",
                                   delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")


crc_belgian_counts <- read_delim("input/crc_data/GSE144735_processed_head.txt",
                                delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")



# Load and check metadata
crc_korean_meta <- read_delim("input/crc_data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",
                              delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")

crc_protocols_meta <- read_delim("input/crc_data/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt",
                                 delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")

crc_belgian_meta <- read_delim("input/crc_data/GSE144735_processed_KUL3_CRC_10X_annotation.txt",
                               delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")





# Create Seurat
seurat <- CreateSeuratObject(counts = Seurat::as.sparse(crc_korean_counts),
                             project = "10X_CRC")
seurat <- Seurat::AddMetaData(seurat, crc_korean_meta)

seurat@meta.data




# Create Seurat
seurat <- CreateSeuratObject(counts = Seurat::as.sparse(crc_protocols_counts),
                             project = "10X_CRC")
seurat <- Seurat::AddMetaData(seurat, crc_protocols_meta)


saveRDS(seurat, "input/crc_data/protocols_seurat.rds")
seurat <- readRDS("input/crc_data/protocols_seurat.rds")


seurat
Seurat::Idents(seurat) <- seurat@meta.data$Cell_type
Idents(seurat)


# SeuratDisk::SaveH5Seurat(seurat, "input/crc_data/crc_korean.h5Seurat")
#
#
# xd <- SeuratDisk::LoadH5Seurat("input/crc_data/crc_korean.h5Seurat")
#
#
#
# unique(seurat@meta.data$Cell_type)
# unique(seurat@meta.data$Cell_subtype)
