

# single cell breast cancer tissue sc data
# use to label spots
# data and markers taken from
# https://www.nature.com/articles/s41523-020-00175-8#Abs1 SI
bc_counts <- read.table("input/sc_bc/GSE118389_norm_data.txt")

bc_seurat <- Seurat::CreateSeuratObject(counts = bc_counts)
FeatureScatter(bc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

bc_seurat <- FindVariableFeatures(bc_seurat, selection.method = "vst", nfeatures = 2000)

# scale
all.genes <- rownames(bc_seurat)
bc_seurat <- ScaleData(bc_seurat, features = all.genes)

# pca
bc_seurat <- RunPCA(bc_seurat, features = VariableFeatures(object = bc_seurat))
ElbowPlot(bc_seurat, ndims = 50)

# cluster
bc_seurat <- FindNeighbors(bc_seurat, dims = 1:20) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:20)

DimPlot(bc_seurat, reduction = "umap")
bc_seurat
# saveRDS(bc_seurat, "input/sc_bc/bc_seurat.rds")


markers <- FindAllMarkers(bc_seurat, logfc.threshold = 1, test.use = "wilcox")
# saveRDS(markers, "input/sc_bc/markers.rds")

top10_gene <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = desc(p_val))
DoHeatmap(bc_seurat, features = top10_gene$gene) + NoLegend()

# General epithelial markers
epith_markers <- c("EPCAM", "EGFR", "CDH1")
FeaturePlot(bc_seurat, features = epith_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = epith_markers)
markers %>% filter(gene %in% epith_markers) %>% filter(p_val_adj < 0.05)


# basal epithelial markers
basal_epith_markerk <-  c("KRT14", "ITGA6", "KRT5", "TP63", "KRT17", "MME")
FeaturePlot(bc_seurat, features = basal_epith_markerk, label = TRUE)
VlnPlot(object = bc_seurat, features = basal_epith_markerk)
markers %>% filter(gene %in% basal_epith_markerk) %>% filter(p_val_adj < 0.05)

# luminal epithelial markers
lum_epith_markers <-  c("KRT8", "KRT18", "KRT19", "FOXA1", "GATA3", "MUC1", "CD24")
FeaturePlot(bc_seurat, features = lum_epith_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = lum_epith_markers)
markers %>% filter(gene %in% lum_epith_markers) %>% filter(p_val_adj < 0.05)
#

# luminal progenitor markers
lum_prog_markers <- c("KIT", "GABRP")
FeaturePlot(bc_seurat, features = lum_prog_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = lum_prog_markers)
markers %>% filter(gene %in% lum_prog_markers) %>% filter(p_val_adj < 0.05)

# stroma markers
stroma_markers <- c("FAP", "COL1A1", "COL3A1", "COL5A1", "ACTA2", "TAGLN",
                    "LUM", "FBLN1", "COL6A3", "COL1A2", "COL6A1", "COL6A2")
FeaturePlot(bc_seurat, features = stroma_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = stroma_markers)
markers %>% filter(gene %in% stroma_markers) %>% filter(p_val_adj < 0.05)
# 5 stroma
# 6 "the cells that express both epithelial and stroma markers are considered
# to indicate the epithelial-mesenchymal transition (EMT) state and are assigned to the epithelial class"
#

# endothelial marker
endo_markers <- c("PECAM1", "VWF", "CDH5", "SELE")
FeaturePlot(bc_seurat, features = endo_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = endo_markers)

# General immune cell marker
FeaturePlot(bc_seurat, features = c("PTPRC"), label = TRUE)
VlnPlot(object = bc_seurat, features = "PTPRC")
# 3,8, 10


# T cell markers
t_markers <- c("CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B")
FeaturePlot(bc_seurat, features = t_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = t_markers)
# clust 3 = T

# B cell markers
b_markers <- c("MS4A1", "CD79A", "CD79B", "BLNK")
FeaturePlot(bc_seurat, features = b_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = b_markers)
# clust 10 = B

# Macrophage markers
macro_markers <- c("CD14", "CD68", "CD163", "CSF1R")
FeaturePlot(bc_seurat, features = macro_markers, label = TRUE)
VlnPlot(object = bc_seurat, features = macro_markers)
# clust 8 = macro

