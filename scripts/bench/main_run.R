# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

sapply(list.files("scripts/pipes/", pattern = ".R", full.names = TRUE), source)

# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat323.rds")
# Fix for NATMI
clust.anns <- c("c0", "c1","c2",
                "c3","c4","c5",
                "c6", "c7", "c8",
                "c9", "c10", "c11", "c12")
names(clust.anns) <- levels(breast_cancer)
breast_cancer <- RenameIdents(breast_cancer, clust.anns)

# Get Omni Resrouces
# source("scripts/utils/get_omnipath.R")
# omni_resources <- get_omni_resources()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <-readRDS("input/omni_resources.rds")

# Get Random DB
source("scripts/utils/shuffle_omnipath.R")
op_random <- shuffle_omnipath(omni_resources$OmniPath)


# Squidpy (benched)

# Append shuffled omni db to omni_resources
omni_resources <- append(list("Random" = op_random,
                              "Default" = NULL),
                         omni_resources)



# SCA
sca_res <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                     seurat_object = breast_cancer,
                     assay = 'SCT',
                     .format = TRUE,
                     s.score = 0,
                     logFC = 0.25
                 ))
saveRDS(sca_res, "output/benchmark/main_run/sca_res.rds")






# CellChat
cellchat_results <- omni_resources %>%
    map(function(x) call_cellchat(x,
                                  seurat_object,
                                  nboot = 100,
                                  exclude_anns = c(),
                                  thresh = 1,
                                  assay = "RNA")) %>%
    setNames(names(omni_resources))


cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = seurat_object,
                                  exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                                  nboot = 100,
                                  thresh = 1,
                                  assay = "RNA")

cellchat_results <- append(cellchat_results,
                           list(cellchat_def = cellchat_default))

# Try with Shuffled DB
cellchat_random <- call_cellchat(op_random,
                                 seurat_object,
                                 nboot = 100,
                                 exclude_anns = c(),
                                 thresh = 1,
                                 assay = "RNA")
