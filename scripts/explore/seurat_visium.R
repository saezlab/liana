library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# InstallData("stxBrain")
# options(timeout=10000)

# 10x Visium
brain <- LoadData("stxBrain", type = "anterior1")

# The visium data from 10x consists of the following data types:
# A spot by gene expression matrix
# An image of the tissue slice (obtained from H&E staining during data acquisition)
# Scaling factors that relate the original high resolution image to the lower resolution image used here for visualization.

# In the Seurat object, the spot by gene expression matrix is similar to a typical "RNA" Assay
# but contains spot level, not single-cell level data.
# The image itself is stored in a new images slot in the Seurat object.
# The images slot also stores the information necessary to associate spots with their physical position on the tissue image.

## Data preprocessing
# The initial preprocessing steps that we perform on the spot by gene expression data are similar to a typical scRNA-seq experiment.
# We first need to normalize the data in order to account for variance in sequencing depth across data points.
# We note that the variance in molecular counts / spot can be substantial for spatial datasets,
# particularly if there are differences in cell density across the tissue. We see substantial heterogeneity here, which requires effective normalization.
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
# These plots demonstrate that the variance in molecular counts across spots is not just technical in nature, but also is dependent on the tissue anatomy


# sctransform builds regularized negative binomial models of gene expressionto account for technical artifacts while preserving biological variance.
# sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)


# Gene expression visualization
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))


# Dimensionality reduction, clustering, and visualization
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)


p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2


## Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue.
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)



# Subset out anatomical regions
# Maybe it would be best to just focus on a specific region of the brain
# rather than the whole brain
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 |
image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2


# saveRDS(cortex, "input/cortex_seruat.rds")
# saveRDS(brain, "input/brain_seruat.rds")
# rm(brain)
# gc()



# load cortex data
cortex <- readRDS("input/cortex_seruat.rds")

# Integration with single-cell data
# At ~50um, spots from the visium assay will encompass the expression profiles of multiple cells
# As such, we deconvolute' each of the spatial voxels to predict the underlying composition of cell types
# we use a single cell reference of 14k cells to transfer cell labels
allen_reference <- readRDS("input/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells
# this speeds up SCTransform dramatically with no loss in performance
library(dplyr)

allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE,
                               conserve.memory = TRUE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# saveRDS(allen_reference, "input/allen_cortex_norm.rds")


# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# saveRDS(cortex, "input/cortex_renormalized.rds")

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)


anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
# saveRDS(anchors, "input/cortex_anchors")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]])
cortex[["predictions"]] <- predictions.assay

# Now we get prediction scores for each spot for each class.
# Of particular interest in the frontal cortex region are the laminar excitatory neurons.
# Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:

DefaultAssay(cortex) <- "SCT"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)



DimPlot(cortex)

# saveRDS(cortex, "input/cortex_deconvoluted.rds")


# Based on these prediction scores, we can also predict cell types whose location is spatially restricted.
# We use the same methods based on marked point processes to define spatially variable features,
# but use the cell type prediction scores as the "marks" rather than gene expression.
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
SpatiallyVariableFeatures(cortex)

# Finally, we show that our integrative procedure is capable of recovering
# the known spatial localization patterns of both neuronal and non-neuronal subsets,
# including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

# saveRDS(cortex, "input/cortex_final.rds")


# Assign cluster names based on markers
DimPlot(cortex)

cortex.anch <- FindAllMarkers(cortex)


FeaturePlot(cortex, features = c("Gfap",
                                 "Cnp"))



# Add spatial info to meta and assay using label transfer ----------------------

# Add spatial to meta
allen_reference <- readRDS("input/allen_cortex_norm.rds")
cortex <- readRDS("input/cortex_renormalized.rds")


anchors <- FindTransferAnchors(reference = allen_reference,
                               query = cortex,
                               normalization.method = "SCT")

labels.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass,
                                  weight.reduction = cortex[["pca"]])
cortex <- AddMetaData(cortex, metadata = labels.assay)

cortex <- SetIdent(cortex, value = cortex@meta.data$predicted.id)
DimPlot(cortex)


# Add spetial to assay
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]])
cortex[["predictions"]] <- predictions.assay



SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

# predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
#                                   weight.reduction = cortex[["pca"]])

cortex@images$anterior1@coordinates

saveRDS(cortex, "input/cortex_final.rds")



