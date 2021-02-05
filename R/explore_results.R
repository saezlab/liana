# Explore Pipe results

# Load prerequisites
library(tidyverse)
library(Seurat)
# library(SeuratData

# load seurat object
seurat_object <- readRDS("input/pbmc3k_processed.rds")
table(Idents(seurat_object))

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


# 2.Connectome
connectome_results <- omni_resources %>%
    map(function(op_resource){
        conn <- call_connectome(op_resource,
                                seurat_object,
                                # optional args passed to createConnectom
                                LR.database = 'custom',
                                min.cells.per.ident = 75,
                                p.values = FALSE,
                                calculate.DOR = FALSE)
    })  %>%
    setNames(names(omni_resources)) %>%
    map(function(conn_res) conn_res %>%
            FilterConnectome(.,
                             min.pct = 0.1,
                             min.z = 0.25,
                             remove.na = TRUE))


connectome_default <- call_connectome_default(seurat_object = seurat_object,
                                              min.cells.per.ident = 75,
                                              p.values = FALSE,
                                              calculate.DOR = FALSE) %>%
    FilterConnectome(.,
                     min.pct = 0.1,
                     min.z = 0.25,
                     remove.na = TRUE)



# 3. NATMI ---------------------------------------------------------------------
