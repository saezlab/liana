SeuratData::AvailableData()
SeuratData::InstallData("cbmc")



# omni_resources <- compile_ligrec(lr_pipeline = TRUE)
omni_resources <- readRDS("input/omni_resources.rds")
op_resources <- list("CellChatDB" = omni_resources$CellChatDB)
op_resources$CellChatDB %<>% filter(str_detect(target_genesymbol, "_"))



# panc8 (cluster)
panc8data <- SeuratData::LoadData("panc8") %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()
summary(as.factor(panc8data@meta.data$celltype))
nlevels(as.factor(panc8data@meta.data$celltype))

panc8data@meta.data %<>%
    mutate(celltype = as.factor(celltype))
panc8data <- subset(panc8data, cells = rownames(panc8data@meta.data))
panc8data <- SetIdent(panc8data, value = panc8data@meta.data$celltype)





# CellChat
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = panc8data,
                        nboot = 10,
                        exclude_anns = c(),
                        thresh = 1,
                        assay = "RNA",
                        .normalize = FALSE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

# Connectome
conn_res <- call_connectome(seurat_object = panc8data,
                            .spatial = FALSE,
                            op_resource = omni_resources$CellPhoneDB,
                            min.cells.per.ident = 1,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            assay = 'RNA',
                            .format = TRUE)

# iTALK
italk_res <- call_italk(op_resource = omni_resources$CellPhoneDB,
                        seurat_object = panc8data,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)


# SCA
sca_res <- call_sca(op_resource = omni_resources$CellPhoneDB,
                    seurat_object = panc8data,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5)
)

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = panc8data,
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             omni_resources = op_resources,
                             .ident = "celltype")







# Local
cbmcdata <- SeuratData::LoadData("cbmc") %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures()

cbmcdata@meta.data %<>%
    mutate(celltype = as.factor(rna_annotations)) %>%
    # remove + as it breaks Squidpy
    mutate(rna_annotations = str_replace_all(rna_annotations, "[+]", "")) %>%
    mutate(rna_annotations = str_replace(rna_annotations, "/", " ")) %>%
    mutate(rna_annotations = as.factor(rna_annotations))


cbmcdata <- subset(cbmcdata, cells = rownames(cbmcdata@meta.data))
cbmcdata <- SetIdent(cbmcdata, value = cbmcdata@meta.data$rna_annotations)

summary(as.factor(cbmcdata@meta.data$rna_annotations))
nlevels(as.factor(cbmcdata@meta.data$rna_annotations))



# CellChat
cc_res <- call_cellchat(op_resource = NULL,
                        seurat_object = cbmcdata,
                        nboot = 10,
                        exclude_anns = c(),
                        thresh = 1,
                        assay = "RNA",
                        .normalize = FALSE,
                        .do_parallel = FALSE,
                        .raw_use = TRUE)

# Connectome
conn_res <- call_connectome(seurat_object = cbmcdata,
                            .spatial = FALSE,
                            op_resource = omni_resources$CellPhoneDB,
                            min.cells.per.ident = 1,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            assay = 'RNA',
                            .format = TRUE)

# iTALK
italk_res <- call_italk(op_resource = omni_resources$CellPhoneDB,
                        seurat_object = cbmcdata,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)


# SCA
sca_res <- call_sca(op_resource = omni_resources$CellPhoneDB,
                    seurat_object = cbmcdata,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5)
)

# Squidpy
squidpy_res <- call_squidpyR(seurat_object = cbmcdata,
                             python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                             omni_resources = op_resources,
                             .ident = "rna_annotations")

