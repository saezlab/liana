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


# 1. CellChat ---------------------------------------------------------------
# Call cellchat and iterate over omni resources
source("scripts/pipes/cellchat_pipe.R")
cellchat_results <- omni_resources %>%
    map(function(x) call_cellchat(x,
                                  seurat_object,
                                  nboot = 1000,
                                  thresh = 0.05)) %>%
    setNames(names(omni_resources))


cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = seurat_object,
                                  exclude_anns = c("ECM-Receptor",
                                                   "Cell-Cell Contact"),
                                  nboot = 1000,
                                  thresh = 0.05)

cellchat_results <- append(cellchat_results,
                           list(cellchat_def = cellchat_default))

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
                                # optional args passed to createConnectom
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

# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources,
#               omni_path = "input/omnipath_NATMI")


# call NATMI
source("scripts/pipes/NATMI_pipe.R")
py_set_seed(1004)
natmi_results <- call_natmi(omni_resources,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/test_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/test_metadata.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/NATMI_test")


# 4. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")
squidpy_results <- call_squidpyR(seurat_object = seurat_object,
                               omni_resources = omni_resources,
                               python_path = "/home/dbdimitrov/anaconda3/bin/python")





# 5. Combine results and plot -------------------------------------------------
# saveRDS(cellchat_results, "output/pbmc3k/cellchat_results.rds")
# saveRDS(connectome_results, "output/pbmc3k/connectome_results.rds")
# saveRDS(natmi_results, "output/pbmc3k/natmi_results.rds")
# saveRDS(squidpy_results, "output/pbmc3k/squidpy_results.rds")

omnipath_LRs <- list("cellchat" = cellchat_results$OmniPath,
                     "connectome" = connectome_results$OmniPath,
                     "natmi" = natmi_results$OmniPath,
                     "squidpy" = squidpy_results$OmniPath) %>%
    map(function(x) x %>% filter(!(.$source==.$target))) # filter autocrine



# Already filtered
omnipath_LRs$cellchat

# Already filtered
omnipath_LRs$connectome



# I would've used edge_specificity, but this is more similar to what they do
# to compare with CellPhoneDB
# they used an arbitrry threshold to filter
omnipath_LRs$natmi <- omnipath_LRs$natmi %>%
    mutate(top5p = ntile(edge_avg_expr, 100)) %>%
    filter(top5p > 95)


omnipath_LRs$squidpy <- omnipath_LRs$squidpy %>%
    filter(pvalue < 0.05)



# Default LRs
default_LRs <- list("cellchat_def" = cellchat_results$cellchat_def,
                    "cellchat_omni" = cellchat_results$CellChatDB,
                    "connectome_def" = connectome_results$connectome_def,
                    "connectome_omni" = connectome_results$connectomeDB2020,
                    "natmi" = natmi_results$lrc2p,
                    "squidpy" = squidpy_results$CellPhoneDB) %>%
    map(function(x) x %>% filter(!(.$source==.$target)))



default_LRs$natmi <- default_LRs$natmi %>%
    mutate(top5p = ntile(edge_avg_expr, 100)) %>%
    filter(top5p > 95)


default_LRs$squidpy <- default_LRs$squidpy %>%
    filter(pvalue < 0.05)


#
library(UpSetR)


