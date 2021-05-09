require(intercell)
require(tibble)
require(purrr)

# Load Data
omni_resources <- readRDS("~/Repos/ligrec_decoupleR/input/omni_resources.rds")
# crc_data <- readRDS("~/Repos/ligrec_decoupleR/input/crc_data/crc_korean_form.rds")


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
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = readRDS("~/Repos/ligrec_decoupleR/input/crc_data/crc_korean_form.rds"),
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = FALSE,
                                   .do_parallel = TRUE,
                                   .raw_use = TRUE
    )) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/crc_res/cellchat_results_loc.rds")
cc_res <- readRDS("~/Repos/ligrec_decoupleR/output/cellchat_local.rds")



cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = readRDS("~/Repos/ligrec_decoupleR/input/crc_data/crc_korean_form.rds"),
                        nboot = 1000,
                        exclude_anns = c(),
                        thresh = 1,
                        assay = "RNA",
                        .normalize = FALSE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       # "pval"=FALSE,
                                       "prob"=TRUE
                                   )))
