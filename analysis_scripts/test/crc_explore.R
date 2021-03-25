# Check Belgian
crc_belgian <- readRDS("input/crc_data/crc_belgian.rds")
crc_korean <- readRDS("input/crc_data/crc_korean.rds")


belgian_clusters <- crc_belgian@meta.data %>%
    unite(Cell_type, Cell_subtype, col = "cc", remove = FALSE, sep=">") %>%
    # Filter Unknown and Unspecified
    filter(!(str_detect("Unknown", Cell_subtype))) %>%
    filter(!(str_starts("Unspecified Plasma", Cell_subtype))) %>%
    # Myoloids_SPP1+A/B into Myoloids_SPP1+ (+ fix typo)
    mutate(cc = if_else(str_detect(cc, "SPP1"),
                        "Myoloids>SPP1+",
                        cc)
           ) %>%
    # Group Healthy Epithelial Cells
    mutate(cc = if_else((str_detect(cc, "Epithelial") & !str_detect(cc, "CMS")),
                        "Epithelial>Healthy",
                        cc
                        )
    ) %>%
    group_by(cc) %>%
    summarise(n())

korean_clusters <- crc_korean@meta.data %>%
    # Filter Unknown Subtypes
    filter(Cell_subtype != "Unknown") %>%
    unite(Cell_type, Cell_subtype, col = "cc", remove = FALSE, sep=">") %>%
    # Myoloids_SPP1+A/B into Myoloids_SPP1+ (+ fix typo)
    mutate(cc = if_else(str_detect(cc, "SPP1"),
                        "Myoloids>SPP1+",
                        cc)
    ) %>%
    # Group Healthy Epithelial Cells
    mutate(cc = if_else((str_detect(cc, "Epithelial") & !str_detect(cc, "CMS")),
                        "Epithelial>Healthy",
                        cc
                        )
           ) %>%
    group_by(cc) %>%
    summarise(n())


setdiff(belgian_clusters$cc, korean_clusters$cc)
setdiff(korean_clusters$cc, belgian_clusters$cc)
# Myeloids_SPP1+A and Myeloids_SPP1+B -> Myeloids_SPP1+
# Epithelial Cells -> Tumour, Stemlike, and Normal
#

xd <- belgian_clusters %>%
    tidyr::separate(cc, into = c("Cluster_types", "Cluster_subtypes"), sep = ">")














