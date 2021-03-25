require(tidyverse)
require(CellChat)
require(Seurat)

source("~/Repos/ligrec_decoupleR/R/cellchat_pipe.R")

# Load Data
omni_resources <- readRDS("~/Repos/ligrec_decoupleR/input/omni_resources.rds")
breast_cancer <- readRDS("~/Repos/ligrec_decoupleR/input/sc_bc/breast_cancer_seurat323.rds")


cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = breast_cancer,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "SCT",
                                   .normalize = FALSE,
                                   .do_parallel = TRUE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "~/Repos/ligrec_decoupleR/output/benchmark/main_run/cellchat_full.rds")


