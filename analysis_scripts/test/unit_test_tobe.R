# omni_resources <- compile_ligrec(lr_pipeline = TRUE)
omni_resources <- readRDS("input/omni_resources.rds")
op_resources <- list("CellChatDB" = omni_resources$CellChatDB)
op_resources$CellChatDB %<>% filter(str_detect(target_genesymbol, "_"))

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

seurat_object <- readRDS("input/testdata.rds")

# CellChat
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = seurat_object,
                        nboot = 10,
                        exclude_anns = c(),
                        thresh = 1,
                        assay = "RNA",
                        .normalize = TRUE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

# Connectome
conn_res <- call_connectome(seurat_object = seurat_object,
                            .spatial = FALSE,
                            op_resource = omni_resources$CellPhoneDB,
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


# SCA
sca_res <- call_sca(op_resource = omni_resources$CellPhoneDB,
                    seurat_object = seurat_object,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5)
)

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = seurat_object,
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             omni_resources = op_resources,
                             n_perms=100,
                             threshold=0.01,
                             seed=as.integer(1004),
                             cluster_key="seurat_annotations")



# NATMI
natmi_results <- call_natmi(omni_resources = op_resources,
                            seurat_object = seurat_object,
                            omnidbs_path = "input/omnipath_NATMI",
                            natmi_path = "NATMI/",
                            em_path = "input/test_em.csv",
                            ann_path = "input/test_metadata.csv",
                            output_path = "output/NATMI_test",
                            .write_data = TRUE,
                            .assay = "RNA"
                            )
xd <- FormatNatmi(output_path,
            names(op_resources),
            TRUE)
