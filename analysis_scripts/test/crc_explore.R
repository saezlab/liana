# Check Belgian
crc_belgian <- readRDS("input/crc_data/crc_belgian.rds")
crc_korean <- readRDS("input/crc_data/crc_korean.rds")

belgian_clusters <- format_crc_meta(crc_belgian) %>%
    pluck("meta.data") %>%
    group_by(Cell_clusters) %>%
    summarise(n())

korean_clusters <- format_crc_meta(crc_korean) %>%
    pluck("meta.data") %>%
    group_by(Cell_subtype) %>%
    summarise(n())




setdiff(belgian_clusters$Cell_clusters, korean_clusters$Cell_clusters)
setdiff(korean_clusters$Cell_clusters, belgian_clusters$Cell_clusters)
# Myeloids_SPP1+A and Myeloids_SPP1+B -> Myeloids_SPP1+
# Epithelial Cells -> Tumour, Stemlike, and Normal
