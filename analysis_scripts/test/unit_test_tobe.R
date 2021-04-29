omni_resources <- compile_ligrec(lr_pipeline = TRUE)
testdata <- SeuratData::LoadData("pbmc3k")
testdata@meta.data <- testdata@meta.data %>%
    filter(seurat_annotations %in% c("NK", "CD8 T", "B")) %>%
    rownames_to_column("id") %>%
    group_by(seurat_annotations) %>%
    top_n(11, nCount_RNA) %>%
    as.data.frame() %>%
    column_to_rownames("id")

testdata <- subset(testdata, cells = rownames(testdata@meta.data))
testdata <- SetIdent(testdata, value = testdata@meta.data$seurat_annotations)


# CellChat
cc_res <- call_cellchat(op_resource = omni_resources$CellPhoneDB,
                    seurat_object = testdata,
                    nboot = 100,
                    exclude_anns = c(),
                    thresh = 0.05,
                    assay = "RNA",
                    .normalize = FALSE,
                    .do_parallel = FALSE,
                    .raw_use = TRUE)

# Connectome (should I filter norm_weight also by p-value?)
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





# NATMI (make it so that the full path is not needed?)

