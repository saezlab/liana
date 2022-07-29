library(SeuratDisk)
library(HDF5Array)
library(scuttle)
library(tidyverse)
library(liana)

h5 <- SeuratDisk::Connect("/media/dbdimitrov/SSDDimitrov/Repos/biopsies/seuratobject.h5Seurat")

h5[["meta.data"]]
print(h5$index())

# Read Sobj
sobj <- SeuratDisk::LoadH5Seurat("/media/dbdimitrov/SSDDimitrov/Repos/biopsies/seuratobject.h5Seurat",
                                 assays = c(SCT = "counts"), #### CHANGE TO RNA!!!!
                                 graphs = FALSE)
sce <- Seurat::as.SingleCellExperiment(sobj)
sce@assays@data@listData$logcounts <- NULL
rm(sobj)

# Basic Feature Filtering
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Save as HDF5 Experiment:
HDF5Array::saveHDF5SummarizedExperiment(sce, "/media/dbdimitrov/SSDDimitrov/Repos/biopsies/hdf5arr/", replace=TRUE)

# Read HDF5
sce <- HDF5Array::loadHDF5SummarizedExperiment("/media/dbdimitrov/SSDDimitrov/Repos/biopsies/hdf5arr/")

# Normalize RNA assay
sce <- scuttle::logNormCounts(sce)

gc()

### Internal Functions for LIANA to filter cell types across samples (For SC)

# Cell types per sample /w at least X cells
# Cell type in at least Y samples
# Proportion of at least Z% of the samples
# sample_col = "Sample"
# condition_col = "Group"
# idents_col = "predicted.annotation.l2"
# min_cells = 20
# min_samples = 3
# min_prop = 0.15


# 1. Filter Samples - e.g. 3 z-scores < of SUM total counts

#' Helper function to filter outlier samples in terms of total counts
#'
#' @param sce SingleCellExperiment
#' @param sample_col Name of the column with sample ids
filter_samples <- function(sce, sample_col){

    pb <- scuttle::aggregateAcrossCells(sce, ids=sce[[sample_col]])

    sample_dist <- scale(colSums(pb@assays@data$counts)) %>%
        as_tibble(rownames="sample",
                  .name_repair = "universal") %>%
        dplyr::rename(zscore = `...1`) %>%
        mutate(keep_total = zscore > -3)

    sample_dist %>%
        ggplot(aes(x = zscore, fill=keep_total)) +
        geom_histogram() +
        scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'red'),
                          guide = "none") +
        geom_vline(aes(xintercept=-3),
                   color="red", linetype="dashed") +
        theme_bw()
    print(sample_dist)

    keep_total <- sample_dist %>% select(sample, keep_total) %>% deframe()
    keep_total <- colData(sce) %>%
        as_tibble(rownames="barcode") %>%
        filter(.data[[sample_col]] %in% names(keep_total[keep_total])) %>%
        pull(barcode)
    # Keep only samples with sum counts within 3 z-scores
    sce <- sce[, keep_total]

    return(sce)
}

sce <- filter_samples(sce, sample_col = "Sample")


# 2. Filter Cell types by min.cell num + min.cell by sample
ctqc <- get_abundance_summary(sce,
                              sample_col = "Sample",
                              idents_col = "predicted.annotation.l2")

# Before filt
plot_abundance_summary(ctqc, ncol=4)


# filter
sce <- filter_nonabundant_celltypes(sce,
                                    ctqc,
                                    sample_col = "Sample",
                                    idents_col = "predicted.annotation.l2"
)


# after filt
ctqc.after <- get_abundance_summary(sce,
                                    sample_col = "Sample",
                                    idents_col = "predicted.annotation.l2")
plot_abundance_summary(ctqc.after, ncol=3)


# + 3a. Filter genes/LR /w edgeR fun? for pseudobulk
# ^liana_pseudobulk



# + 3a. Filter Proportions of Cell types and LR Expression samples (Untargetted)
# ^MOFAphone


# --- Test Keep only IgA and TRL
sce$Group %>% unique()

groups_of_interest <- colData(sce) %>%
    as_tibble(rownames="barcodes") %>%
    filter(Group %in% c("IgA", "CTRL")) %>%
    pull("barcodes")
sce <- sce[,groups_of_interest]


# Run LIANA by Sample ----
context_df_dict <- liana_bysample(sce = sce,
                                  sample_col = "Sample",
                                  condition_col = "Group",
                                  idents_col = "predicted.annotation.l2",
                                  method="sca")
saveRDS(context_df_dict, "/media/dbdimitrov/SSDDimitrov/Repos/biopsies/context_df_dict.RDS")

