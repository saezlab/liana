library(Seurat)
library(tidyverse)
# library("hdf5r")
library(SeuratDisk)
library(ggplot2)
library(patchwork)


# remotes::install_github("mojaveazure/seurat-disk")
breast_cancer <- Load10X_Spatial("input/human_breast_cancer/",
                filename = "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5")


# install.packages("Seurat")
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


?SCTransform

breast_cancer
# saveRDS(breast_cancer, "input/sc_bc/breast_cancer_seurat.rds")
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat.rds")

GetAssay(breast_cancer)[[]]


# Annotate ----
markers <- FindAllMarkers(breast_cancer, logfc.threshold = 0.25, test.use = "wilcox")



# General epithelial markers
epith_markers <- c("EPCAM", "EGFR", "CDH1")
FeaturePlot(breast_cancer, features = epith_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = epith_markers)
markers %>% filter(gene %in% epith_markers) %>% filter(p_val_adj < 0.05)


# basal epithelial markers
basal_epith_markerk <-  c("KRT14", "ITGA6", "KRT5", "TP63", "KRT17", "MME")
FeaturePlot(breast_cancer, features = basal_epith_markerk, label = TRUE)
VlnPlot(object = breast_cancer, features = basal_epith_markerk)
markers %>% filter(gene %in% basal_epith_markerk) %>% filter(p_val_adj < 0.05)

# luminal epithelial markers
lum_epith_markers <-  c("KRT8", "KRT18", "KRT19", "FOXA1", "GATA3", "MUC1", "CD24")
FeaturePlot(breast_cancer, features = lum_epith_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = lum_epith_markers)
markers %>% filter(gene %in% lum_epith_markers) %>% filter(p_val_adj < 0.05)
#

# luminal progenitor markers
lum_prog_markers <- c("KIT", "GABRP")
FeaturePlot(breast_cancer, features = lum_prog_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = lum_prog_markers)
markers %>% filter(gene %in% lum_prog_markers) %>% filter(p_val_adj < 0.05)

# stroma markers
stroma_markers <- c("FAP", "COL1A1", "COL3A1", "COL5A1", "ACTA2", "TAGLN",
                    "LUM", "FBLN1", "COL6A3", "COL1A2", "COL6A1", "COL6A2")
FeaturePlot(breast_cancer, features = stroma_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = stroma_markers)
markers %>% filter(gene %in% stroma_markers) %>% filter(p_val_adj < 0.05)

# endothelial marker
endo_markers <- c("PECAM1", "VWF", "CDH5", "SELE")
FeaturePlot(breast_cancer, features = endo_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = endo_markers)

# General immune cell marker
FeaturePlot(breast_cancer, features = c("PTPRC"), label = TRUE)
VlnPlot(object = breast_cancer, features = "PTPRC")


# T cell markers
t_markers <- c("CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B")
FeaturePlot(breast_cancer, features = t_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = t_markers)

# B cell markers
b_markers <- c("MS4A1", "CD79A", "CD79B", "BLNK")
FeaturePlot(breast_cancer, features = b_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = b_markers)

# Macrophage markers
macro_markers <- c("CD14", "CD68", "CD163", "CSF1R")
FeaturePlot(breast_cancer, features = macro_markers, label = TRUE)
VlnPlot(object = breast_cancer, features = macro_markers)


# Annotate unknown clusters with scCATCH
# devtools::install_github('ZJUFanLab/scCATCH')
breast_cancer@assays$RNA <- breast_cancer@assays$Spatial
library(scCATCH)
clu_markers <- findmarkergenes(breast_cancer,
                               species = "Human",
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = "Breast Cancer", # "Breast Cancer"
                               tissue = "Breast",
                               cell_min_pct = 0.25,
                               logfc = 0.5,
                               pvalue = 0.05)


clu_ann <- scCATCH(clu_markers,
                   species = "Human",
                   cancer = "Breast Cancer",
                   tissue = "Breast")



# immune cell atlas (Seurat object) ------
# This atlas works well for Immune Cell clusters, but we need to extend it
# to also include other cell types
cancer_reference <- readRDS("input/cancer_atlas/TICAtlas_downsampled_1000.rds")

cancer_reference <- SCTransform(cancer_reference, verbose = FALSE,
                                conserve.memory = TRUE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
DimPlot(cancer_reference, group.by = "cell_type", label = TRUE)

# keep only breast cancer cells?
breast_atlas <- cancer_reference@meta.data %>% filter(source == "breast")
breast_atlas <- cancer_reference[ ,rownames(breast_atlas)]

DimPlot(breast_atlas, group.by = "cell_type", label = TRUE)



anchors <- FindTransferAnchors(reference = cancer_reference,
                               query = breast_cancer,
                               normalization.method = "SCT")

# add to meta

labels.assay <- TransferData(anchorset = anchors, refdata = cancer_reference$cell_type,
                             weight.reduction = breast_cancer[["pca"]],
                             dims=1:50)
breast_cancer <- AddMetaData(breast_cancer, metadata = labels.assay)

# Add spetial to assay
predictions.assay <- TransferData(anchorset = anchors,
                                  refdata = cancer_reference$cell_type,
                                  prediction.assay = TRUE,
                                  weight.reduction = breast_cancer[["pca"]],
                                  dims=1:50)
breast_cancer[["predictions"]] <- predictions.assay



breast_cancer@meta.data %>%
    filter(seurat_clusters==6)

