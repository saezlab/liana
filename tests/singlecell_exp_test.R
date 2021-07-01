# library(scRNAseq) # more datasets
require(tidymodels)
require(SingleCellExperiment)
require(scuttle) # util funcs
require(scran)
# require(scater) # visualize

# try on my data
testdata <- readRDS("inst/testdata/input/testdata.rds")

# Convert to SCE
test_sce <- Seurat::as.SingleCellExperiment(testdata)
colLabels(test_sce) <- Seurat::Idents(testdata)


test_summ <- scuttle::summarizeAssayByGroup(test_sce,
                                            ids = colLabels(test_sce))
test_summ@colData
test_summ@assays@data$mean # gene mean across cell types
test_summ@assays@data$prop.detected # gene prop
# test_summ@assays@data$num.detected # num times detected

# (global) Avg Expr by gene
scuttle::calculateAverage(test_sce)

# Pseudo-bulk
pseudobulk <- scuttle::aggregateAcrossCells(test_sce, colLabels(test_sce))


test_t <- scran::findMarkers(test_sce,
                             groups = colLabels(test_sce),
                             direction = "any",
                             full.stats = TRUE,
                             test.type = "wilcox")
test_t@listData





# t-test by cluster (all genes together)
test_summ@assays@data$mean %>%
    as_tibble() %>%
    select("B", "NK") %>%
    pivot_longer(names_to = "celltype",
                 cols = c("B", "NK"),
                 values_to = "counts")  %>%
    t_test(counts ~ celltype,
           order = c("B", "NK"))


