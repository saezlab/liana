#' Helper function to convert CRC data to sparse Seurat
#' @param counts_loc counts location
#' @param meta_loc metadata location
#' @param save_loc location to save sparsified Seurat object
#' @return Return a Sparse Seurat Object
sparsify_to_seurat <- function(counts_loc, meta_loc, save_loc){
    counts <- read_delim(counts_loc,
                         delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    meta <- read_delim(meta_loc,
                       delim = "\t") %>%
        as.data.frame() %>%
        column_to_rownames("Index")

    crc_seurat <- CreateSeuratObject(counts = Seurat::as.sparse(counts),
                       project = "10X_CRC") %>%
        Seurat::AddMetaData(meta) %>%
        Seurat::NormalizeData() %>%
        Seurat::FindVariableFeatures() %>%
        saveRDS(., save_loc)
}


#' Helper function to format metadata for CRC data sets
#' @param crc_seurat CRC Seurat Object
#' @return Return a formatted Metadata for CCC inference
#' @export
format_crc_meta <- function(crc_seurat){
    crc_seurat@meta.data <- crc_seurat@meta.data %>%
        # keep only Tumour samples
        filter(str_detect(orig.ident, "-T")) %>%
        unite(Cell_type, Cell_subtype, col = "Cell_clusters", remove = FALSE, sep="_") %>%
        # Filter Mast, Unknown, and Unspecified
        filter(!(str_detect("Mast cells_Mast cells", Cell_clusters))) %>%
        filter(!(str_detect("Unknown", Cell_subtype))) %>%
        filter(!(str_starts("Unspecified Plasma", Cell_clusters))) %>%
        # Myoloids_SPP1+A/B into Myoloids_SPP1+ (+ fix typo)
        mutate(Cell_clusters = if_else(str_detect(Cell_clusters, "SPP1"),
                                       "Myoloids_SPP1+",
                                       Cell_clusters)
        ) %>%
        # remove + as it breaks Squidpy
        mutate(Cell_subtype = str_replace_all(Cell_subtype, "[+]", "")) %>%
        mutate(Cell_subtype = factor(Cell_subtype))

    crc_seurat <- subset(crc_seurat, cells = rownames(crc_seurat@meta.data))
    crc_seurat <- SetIdent(crc_seurat, value = crc_seurat@meta.data$Cell_subtype)

    return(crc_seurat)
}
