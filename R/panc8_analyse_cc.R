require(intercell)
require(tibble)
require(magrittr)
require(purrr)

setwd("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/")

# Load Data and Format data
# panc8data <- SeuratData::LoadData("panc8") %>%
#     Seurat::NormalizeData() %>%
#     Seurat::FindVariableFeatures()
# panc8data@meta.data %<>%
#     mutate(celltype = as.factor(celltype))
# panc8data <- subset(panc8data, cells = rownames(panc8data@meta.data))
# panc8data <- SetIdent(panc8data, value = panc8data@meta.data$celltype)
# saveRDS(panc8data, "input/panc8_seurat.rds")

panc8data <- readRDS("input/panc8_seurat.rds")

# Get Full Omni Resources
# omni_resources <- compile_ligrec(lr_pipeline = TRUE)
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")


# 6. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = panc8data,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = FALSE,
                                   .do_parallel = FALSE,
                                   .raw_use = TRUE
    )) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/panc8_res/cellchat_results.rds")
