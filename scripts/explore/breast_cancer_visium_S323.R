library(tidyverse)
library(Seurat)

breast_cancer <- Load10X_Spatial("input/human_breast_cancer/",
                                 filename = "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5")


# remotes::install_version("Seurat", version = "3.2.3")
# library(Seurat)


plot1 <- VlnPlot(breast_cancer, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(breast_cancer, features = "nCount_Spatial") + theme(legend.position = "right")
plot1
plot2


breast_cancer <- breast_cancer %>%
    SCTransform(assay = "Spatial", verbose = FALSE) %>%
    RunPCA(assay = "SCT", verbose = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>% # elbow seems to be at ~15
    FindClusters(resolution = 0.8, verbose = TRUE) %>%
    RunUMAP(reduction = "pca", dims = 1:20)


p1 <- DimPlot(breast_cancer, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(breast_cancer, label = TRUE, label.size = 3)
p1 + p2

breast_cancer <- FindSpatiallyVariableFeatures(breast_cancer, assay = "SCT",
                                               features = VariableFeatures(breast_cancer)[1:1000],
                                               selection.method = "markvariogram")
saveRDS(breast_cancer, "input/sc_bc/breast_cancer_seurat323.rds")


