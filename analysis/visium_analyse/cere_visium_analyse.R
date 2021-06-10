# Load data
seurat_object <- Load10X_Spatial("input/cerebellum/",
                                 filename = "Parent_Visium_Human_Cerebellum_filtered_feature_bc_matrix.h5")

# Quality Control Check
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(seurat_object, features = "nCount_Spatial") + theme(legend.position = "right")

# Normalize and Clust
seurat_object <- seurat_object %>%
    SCTransform(assay = "Spatial", verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:15) %>% # elbow seems to be at ~8
    FindClusters(resolution = 0.8, verbose = TRUE) %>%
    RunUMAP(reduction = "pca", dims = 1:15)

# Rename clusters
clust.anns <- str_glue("c{levels(seurat_object)}")
names(clust.anns) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, clust.anns)
seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(seurat_clusters = as.factor(str_glue("c{seurat_clusters}")))
Idents(seurat_object)

# Check clust
DimPlot(seurat_object, reduction = "umap", label = TRUE)
SpatialDimPlot(seurat_object,
               label = TRUE,
               label.size = 3)

# Get Save Object
saveRDS(seurat_object, "input/visium_converted/cere_visium.rds")
