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

# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = panc8data,
                                 omni_resources = omni_resources,
                                 python_path = "/net/data.isilon/ag-saez/bq_ddimitrov/SOFTWARE/miniconda3/envs/ligrec/bin/python3.8",
                                 .ident = "celltype")
saveRDS(squidpy_results, "output/panc8_res/squidpy_results.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = panc8data,
                            wd_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/",
                            omnidbs_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/omni_natmi",
                            natmi_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/NATMI/",
                            em_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/panc8data_counts.csv",
                            ann_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/panc8data_ann.csv",
                            output_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/output/panc8_natmi/",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA",
                            .num_cor = 24
)
saveRDS(natmi_results, "output/panc8_res/natmi_results.rds")


# 3. SCA ----------------------------------------------------------------------
sca_results <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                 panc8data,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5)
        ))
saveRDS(sca_results, "output/panc8_res/sca_results.rds")


# 4. Connectome ----------------------------------------------------------------
conn_results <- omni_resources %>%
    map(function(db)
        call_connectome(seurat_object = panc8data,
                        .spatial = FALSE,
                        op_resource = db,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = TRUE,
                        assay = 'RNA'))
saveRDS(conn_results, "output/panc8_res/conn_results.rds")


# 5. iTALK ---------------------------------------------------------------------
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   seurat_object = panc8data,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/panc8_res/italk_results.rds")



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
