library(data.table)
library(Matrix)

crc_counts <- read_delim("input/crc_data/crc_processed.txt",
                         delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")


# Load and check metadata
crc_korean_meta <- read_delim("input/crc_data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",
                       delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")

korean_clusters <- crc_korean_meta %>%
    group_by(Cell_type, Cell_subtype) %>%
    summarise(n())


crc_protocols_meta <- read_delim("input/crc_data/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt",
                              delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")

protocols_clusters <- crc_protocols_meta %>%
    group_by(Cell_type, Cell_subtype) %>%
    summarise(n())


crc_belgian_meta <- read_delim("input/crc_data/GSE144735_processed_KUL3_CRC_10X_annotation.txt",
                                 delim = "\t") %>%
    as.data.frame() %>%
    column_to_rownames("Index")


belgian_clusters <- crc_belgian_meta %>%
    group_by(Cell_type, Cell_subtype) %>%
    summarise(n())





# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x
                     , mart = human, attributesL = c("mgi_symbol"),
                     martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}

unique(omni_resources$OmniPath)
genesV2[, 2]


genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                 values = omni_resources$Random$source,
                 mart = human, attributesL = c("mgi_symbol"),
                 martL = mouse, uniqueRows=T)

unique(genesV2)
