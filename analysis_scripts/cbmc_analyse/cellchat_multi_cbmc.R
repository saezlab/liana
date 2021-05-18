require(intercell)
require(tibble)
require(purrr)
library(future)

# Load Data
omni_resources <- readRDS("~/Repos/ligrec_decoupleR/input/omni_resources.rds")
cbmcdata <- readRDS("input/cbmc_seurat.rds")


cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = cbmcdata,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = FALSE,
                                   .do_parallel = FALSE,
                                   .raw_use = TRUE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/cbmc_res/cellchat_results.rds")
