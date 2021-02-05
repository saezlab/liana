# Explore Pipe results

# Load prerequisites
library(tidyverse)
library(Seurat)

# load seurat object
seurat_object <- readRDS("input/pbmc3k_processed.rds")


# Get Omni Resrouces
source("./R/utils/get_omnipath.R")
omni_resources <- get_omni_resources()

# 1. CellChat ---------------------------------------------------------------
# Call cellchat and iterate over omni resources
cellchat_results <- omni_resources %>%
    map(function(x) call_cellchat(x,
                                  seurat_object,
                                  thresh = 0.05)) %>%
    setNames(names(omni_resources))


cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = seurat_object,
                                  use_omni_resource = FALSE,
                                  exclude_anns = c("ECM-Receptor",
                                                   "Cell-Cell Contact"),
                                  thresh = 0.05)

# Conclusions:
# Good idea, takes multiple things into considerations that other packages ignore
# cofactors, etc, imputes via PPI, etc.
# However, very difficult to extend Resource, Slow.
# When using OmniPath CellChatDB, we get a slightly different number of
# "significant" hits than when using the original CellChatDB.
# This is liekly due to difference in number of interactions after filtering by annotation
# I opened a Github issue about this and confirmed that what I've done is correct.




