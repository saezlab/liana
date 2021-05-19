# omni_resources <- compile_ligrec(lr_pipeline = TRUE)
omni_resources <- readRDS("input/omni_resources.rds")
op_resources <- list("CellChatDB" = omni_resources$CellChatDB)
op_resources$CellChatDB %<>% filter(str_detect(target_genesymbol, "_"))

testdata <- SeuratData::LoadData("pbmc3k")
testdata@meta.data <- testdata@meta.data %>%
    filter(seurat_annotations %in% c("NK", "CD8 T", "B"))# %>%
    # rownames_to_column("id") %>%
    # group_by(seurat_annotations) %>%
    # top_n(11, nCount_RNA) %>%
    # as.data.frame() %>%
    # column_to_rownames("id")
testdata <- testdata %>% FindVariableFeatures()
testdata <- subset(testdata, cells = rownames(testdata@meta.data))
testdata <- SetIdent(testdata, value = testdata@meta.data$seurat_annotations)


# CellChat
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = testdata,
                        nboot = 100,
                        exclude_anns = c(),
                        thresh = 1,
                        assay = "RNA",
                        .normalize = TRUE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

# Connectome
conn_res <- call_connectome(seurat_object = testdata,
                            .spatial = FALSE,
                            op_resource = omni_resources$CellPhoneDB,
                            min.cells.per.ident = 1,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            assay = 'RNA',
                            .format = TRUE)

# iTALK
italk_res <- call_italk(op_resource = omni_resources$CellPhoneDB,
                        seurat_object = testdata,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)


# SCA
sca_res <- call_sca(op_resource = omni_resources$CellPhoneDB,
                    seurat_object = testdata,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5)
                    )

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = testdata,
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             omni_resources = op_resources,
                             n_perms=10000,
                             threshold=0.1,
                             seed=as.integer(1004),
                             cluster_key="seurat_annotations")


testdata <- SeuratData::LoadData("pbmc3k") %>% FindVariableFeatures()
seurat_object = testdata
omni_resources <- readRDS("input/omni_resources.rds")
op_resources <- list("CellChatDB" = omni_resources$CellChatDB)
op_resources$CellChatDB %<>% filter(str_detect(target_genesymbol, "_"))
python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python"
omni_resources = op_resources
.ident = "seurat_annotations"
.seed = 1004

# NATMI (make it so that the full path is not needed?)


xd <- py$squidpy_results$meta

xd[[1]][[2]]$metadata

# Real data test --------
# Test with real data
# readRDS("input/crc_data/crc_belgian.rds") %>%
#     format_crc_meta() %>%
#     saveRDS("input/crc_data/crc_belgian_form.rds")
# crc_belgian <- readRDS("input/crc_data/crc_belgian_form.rds")

# CellChat
start <- Sys.time()
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = readRDS("input/crc_data/crc_belgian_form.rds"),
                        nboot = 100,
                        exclude_anns = c(),
                        thresh = 0.05,
                        assay = "RNA",
                        .normalize = FALSE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)
end <- Sys.time()

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = readRDS("input/crc_data/crc_belgian_form.rds"),
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             omni_resources = op_resources,
                             .ident = "Cell_subtype")


