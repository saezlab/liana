omni_resources <- compile_ligrec(lr_pipeline = TRUE)
testdata <- SeuratData::LoadData("pbmc3k")
testdata@meta.data <- testdata@meta.data %>%
    filter(seurat_annotations %in% c("NK", "CD8 T", "B")) # %>%
    # rownames_to_column("id") %>%
    # group_by(seurat_annotations) %>%
    # top_n(11, nCount_RNA) %>%
    # as.data.frame() %>%
    # column_to_rownames("id")

testdata <- subset(testdata, cells = rownames(testdata@meta.data))
testdata <- SetIdent(testdata, value = testdata@meta.data$seurat_annotations)

cc <- call_cellchat(op_resource = NULL,
                    seurat_object = testdata,
                    nboot = 100,
                    exclude_anns = c(),
                    thresh = 0.05,
                    assay = "RNA",
                    .normalize = FALSE,
                    .do_parallel = FALSE,
                    .raw_use = TRUE)



seurat_object <- testdata
labels <- Idents(seurat_object)
meta <- data.frame(group = labels, row.names = names(labels))


library(CellChat)

cellchat.omni <- createCellChat(object =
                                   GetAssayData(seurat_object,
                                                assay = "RNA",
                                                slot = "data"),
                           meta = meta,
                           group.by = "group")

cellchat.omni@DB <- CellChatDB.human

cellchat.omni <- subsetData(cellchat.omni)

cellchat.omni <- identifyOverExpressedGenes(cellchat.omni)
cellchat.omni <- identifyOverExpressedInteractions(cellchat.omni)
cellchat.omni <- projectData(cellchat.omni, PPI.human)


cellchat.omni <- computeCommunProb(cellchat.omni, raw.use = TRUE)
cellchat.omni <- filterCommunication(cellchat.omni, min.cells = 1)

df.net <- subsetCommunication(cellchat.omni, thresh = 0.05)

