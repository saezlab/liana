library(Seurat)
library(tidyverse)


# 1. Check if cells close together have a higher expression correlation
cortex <- readRDS("input/cortex_final.rds")

predicted_labels <- cortex@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode, predicted.id)

# check number of spots by label
x <- as.list(summary(as.factor(predicted_labels$predicted.id)))
ok_labels <- names(x[x>10])

# get coordinates with labels and euc dist
# I should look into pairwise euc distance
coords <- cortex@images$anterior1@coordinates %>%
  rownames_to_column("barcode") %>%
  left_join(predicted_labels) %>%
  mutate(euc_dist = map2(row, col, function(x,y){
    sqrt((x - y)^2)
  }))




# We use coordinates, split by cell type


GetAssayData(cortex)

# Get Assay to matrix
assay_matrix <- as.matrix(GetAssayData(cortex))



spots_by_label <- assay_matrix %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  left_join(predicted_labels) # %>%
  # group_by("predicted.id") %>%
  # group_split()



# get most variable features
cortex <- FindVariableFeatures(cortex,
                               nfeatures = 1000)

var_features <- Seurat::VariableFeatures(cortex)

# get average expression of most variable features
avg_expression <- AverageExpression(cortex, assays = c("SCT",
                                                       "Spatial"),
                                    features = var_features
                                    )

# get gene expression correlation between celltypes
clust_corr <-as.data.frame(cor(avg_expression$SCT, method = "pearson"))


# filter non-abundant spots and convert to matrix
clust_corrf <- clust_corr %>% 
  as.data.frame() %>%
  rownames_to_column("celltype") %>%
  select(celltype, ok_labels) %>%
  filter(celltype %in% ok_labels) %>%
  arrange(celltype) %>%
  column_to_rownames("celltype") %>%
  as.matrix()


# plot correlation heatmap
pheatmap::pheatmap(clust_corrf)


# LR - Distance Correlation
