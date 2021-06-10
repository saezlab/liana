# Unit test for omni_resources
omni_resources <- compile_ligrec(lr_pipeline = TRUE)
omni_resources <- readRDS("data/input/omni_resources.rds")
op_resources <- omni_resources[1]

# seurat_object <- SeuratData::LoadData("pbmc3k")
# seurat_object@meta.data <- seurat_object@meta.data %>%
#     filter(seurat_annotations %in% c("NK", "CD8 T", "B")) %>%
#     rownames_to_column("id") %>%
#     group_by(seurat_annotations) %>%
#     top_n(30, nCount_RNA) %>%
#     as.data.frame() %>%
#     column_to_rownames("id")
#
# seurat_object <- Seurat::CreateSeuratObject(counts = seurat_object@assays$RNA@counts,
#                                             meta.data = seurat_object@meta.data)
# seurat_object@meta.data %<>% na.omit()
# seurat_object@meta.data$seurat_annotations %<>% factor()
# seurat_object <- subset(seurat_object, cells = rownames(seurat_object@meta.data))
# seurat_object <- SetIdent(seurat_object, value = seurat_object@meta.data$seurat_annotations)
# seurat_object %<>% FindVariableFeatures()
# saveRDS(seurat_object, "input/testdata.rds")

seurat_object <- readRDS("data/input/testdata.rds")

# Unit tests for each method

# CellChat
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = seurat_object,
                        nboot = 10,
                        exclude_anns = NULL,
                        thresh = 1,
                        assay = "RNA",
                        .normalize = TRUE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

# Connectome
conn_res <- call_connectome(seurat_object = seurat_object,
                            op_resource = omni_resources$CellPhoneDB,
                            .spatial = FALSE,
                            min.cells.per.ident = 1,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            assay = 'RNA',
                            .format = TRUE)

# iTALK
italk_res <- call_italk(op_resource = omni_resources$CellPhoneDB,
                        seurat_object = seurat_object,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)



# NATMI
natmi_results <- call_natmi(op_resource = op_resources,
                            seurat_object = seurat_object,
                            omnidbs_dir = "omnipath_NATMI",
                            expr_file = "em.csv",
                            meta_file = "metadata.csv",
                            output_dir = "NATMI_test",
                            assay = "RNA",
                            num_cor = 4,
                            .format = TRUE,
                            .write_data = TRUE,
                            .seed = 1004,
                            .natmi_path = NULL)


# SCA
sca_res <- call_sca(op_resource = NULL,
                    seurat_object = seurat_object,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5)
)

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = seurat_object,
                             op_resource = op_resources,
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             cluster_key="seurat_annotations",
                             n_perms=1000,
                             threshold=0.01,
                             seed=as.integer(1004))

# check available resources
get_lr_resources()

# LIANA
testrun <- liana_wrap(seurat_object,
                      method = c('italk', 'sca', 'cellchat', 'connectome', 'squidpy', 'natmi'),
                      resource = c('OmniPath', 'connectomeDB2020'),
                      cellchat.defaults = (list(
                              nboot = 10,
                              exclude_anns = NULL,
                              thresh = 1,
                              assay = "RNA",
                              .normalize = FALSE,
                              .do_parallel = FALSE,
                              .raw_use = TRUE
                          )))
