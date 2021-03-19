# Load data
breast_cancer <- Load10X_Spatial("input/human_breast_cancer/",
                                 filename = "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5")

# Quality Control Check
breast_cancer[["percent.mt"]] <- PercentageFeatureSet(breast_cancer, pattern = "^MT-")
VlnPlot(breast_cancer, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(breast_cancer, features = "nCount_Spatial") + theme(legend.position = "right")

breast_cancer <- subset(breast_cancer,
       subset = nCount_Spatial > 500 &
           nFeature_Spatial > 200 &
           percent.mt  < 5)
# Check if improved
VlnPlot(breast_cancer, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(breast_cancer, features = "nCount_Spatial") + theme(legend.position = "right")

# Normalize and Clust
breast_cancer <- breast_cancer %>%
    SCTransform(assay = "Spatial", verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>% # elbow seems to be at ~15
    FindClusters(resolution = 0.8, verbose = TRUE) %>%
    RunUMAP(reduction = "pca", dims = 1:20)


# Rename clusters
clust.anns <- str_glue("c{levels(breast_cancer)}")
names(clust.anns) <- levels(breast_cancer)
breast_cancer <- RenameIdents(breast_cancer, clust.anns)
breast_cancer@meta.data <- breast_cancer@meta.data %>%
    mutate(seurat_clusters = as.factor(str_glue("c{seurat_clusters}")))
Idents(breast_cancer)

# Check clust
DimPlot(breast_cancer, reduction = "umap", label = TRUE)
SpatialDimPlot(breast_cancer,
               label = TRUE,
               label.size = 3)

# Get Save Object
saveRDS(breast_cancer, "input/sc_bc/breast_cancer_seurat323.rds")


# Get NES
# Need to convert NES test to function
# if we want to share the results

