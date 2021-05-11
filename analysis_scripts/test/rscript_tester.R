require(intercell)
require(tidyverse)

crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
    format_crc_meta()

cellchat_results <- call_cellchat(op_resource = NULL,
                                  seurat_object = readRDS("input/crc_data/crc_belgian_form.rds"),
                                  nboot = 100,
                                  # exclude_anns = c("Secreted Signaling"),
                                  thresh = 0.05,
                                  assay = "RNA",
                                  .normalize = FALSE,
                                  .do_parallel = FALSE,
                                  .raw_use = TRUE)
saveRDS(cellchat_results, "~/Repos/ligrec_decoupleR/output/cellchat_local.rds")
cc_res <- readRDS("~/Repos/ligrec_decoupleR/output/cellchat_local.rds")
