#### Explore Pipe results

# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

# load python repl
reticulate::repl_python()

# load test seurat object
seurat_object <- readRDS("input/pbmc3k_processed.rds")
table(Idents(seurat_object))

# Get Omni Resrouces
source("scripts/utils/get_omnipath.R")
omni_resources <- get_omni_resources()

# Get Random DB
source("scripts/utils/shuffle_omnipath.R")
op_random <- shuffle_omnipath(omni_resources$OmniPath)


# 1. CellChat ---------------------------------------------------------------
# Call cellchat and iterate over omni resources
source("scripts/pipes/cellchat_pipe.R")
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

# Conclusions:
# Good idea, takes multiple things into considerations that other packages ignore
# cofactors, etc, imputes via PPI, etc.
# However, very difficult to extend Resource, VERY Slow.
# When using OmniPath CellChatDB, we get a slightly different number of
# "significant" hits than when using the original CellChatDB.
# This is liekly due to difference in number of interactions after filtering by annotation
# I opened a Github issue about this and confirmed that what I've done is correct.


# 2.Connectome
source("scripts/pipes/connectome_pipe.R")
connectome_results <- omni_resources %>%
    map(function(op_resource){
        conn <- call_connectome(op_resource,
                                seurat_object,
                                # optional args passed to createConnectome
                                LR.database = 'custom',
                                min.cells.per.ident = 1,
                                p.values = TRUE,
                                calculate.DOR = FALSE)
    })  %>%
    setNames(names(omni_resources))

connectome_default <- call_connectome(op_resource = NULL,
                                      seurat_object = seurat_object,
                                      min.cells.per.ident = 1,
                                      p.values = TRUE,
                                      calculate.DOR = FALSE,
                                      .format = TRUE)


connectome_results <- append(connectome_results,
                           list(connectome_def = connectome_default))

connectome_random <- call_connectome(op_random,
                                     seurat_object,
                                     LR.database = 'custom',
                                     min.cells.per.ident = 1,
                                     p.values = TRUE,
                                     calculate.DOR = FALSE)




# 3. NATMI ---------------------------------------------------------------------
# Extract data from Seurat Object
# write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object,
#                                             assay = "RNA",
#                                             slot = "data"))) - 1),
#           file = "input/test_em.csv",
#           row.names = TRUE)
# write.csv(Idents(object = seurat_object)  %>%
#               enframe(name="barcode", value="annotation"),
#           file = "input/test_metadata.csv",
#           row.names = FALSE)

source("scripts/pipes/NATMI_pipe.R")

# save OmniPath Resource to NATMI format
omni_plus_random <- append(omni_resources,
                           list("Random" = op_random))
omni_to_NATMI(omni_plus_random,
              omni_path = "input/omnipath_NATMI")


# call NATMI
py_set_seed(1004)
natmi_results <- call_natmi(omni_plus_random,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/test_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/test_metadata.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/NATMI_test")


# 4. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")
squidpy_results <- call_squidpyR(seurat_object = seurat_object,
                               omni_resources = omni_resources,
                               python_path = "/home/dbdimitrov/anaconda3/bin/python",
                               ident = "seurat_clusters")



