require(tidyverse)
require(CellChat)
require(Seurat)
require(logger)

source("~/Repos/ligrec_decoupleR/R/cellchat_pipe.R")

# Load Data
omni_resources <- readRDS("~/Repos/ligrec_decoupleR/input/omni_resources.rds")
crc_data <- readRDS("~/Repos/ligrec_decoupleR/input/crc_data/crc_korean_form.rds")


# cellchat_results <- omni_resources %>%
#     map(function(db) call_cellchat(op_resource = db,
#                                    seurat_object = crc_data,
#                                    nboot = 100,
#                                    exclude_anns = c(),
#                                    thresh = 1,
#                                    assay = "RNA",
#                                    .normalize = TRUE,
#                                    .do_parallel = TRUE)) %>%
#     setNames(names(omni_resources))
cellchat_results <- call_cellchat(op_resource = NULL,
                                   seurat_object = crc_data,
                                   nboot = 100,
                                   exclude_anns = c("Secreted Signaling"),
                                   thresh = 1,
                                   assay = "RNA",
                                  .normalize = FALSE,
                                  .do_parallel = FALSE,
                                  .raw_use = TRUE)
saveRDS(cellchat_results, "~/Repos/ligrec_decoupleR/output/cellchat_local.rds")


cr <- readRDS("~/Repos/ligrec_decoupleR/output/cellchat_local.rds")
