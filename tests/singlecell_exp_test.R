# library(scRNAseq) # more datasets
require(tidyverse)
require(tidymodels)
require(SingleCellExperiment)
require(scuttle) # util funcs
require(scran)
# require(scater) # visualize

sce <- ZeiselBrainData()
colLabels(sce) <- sce@colData$level1class


# Summarize by group
xx <- summarizeAssayByGroup(sce, colLabels(sce))
xx@colData
xx@assays@data$mean # gene mean across cell types
xx@assays@data$prop.detected # gene prop
xx@assays@data$num.detected # num times detected

# try on my data
testdata <- readRDS("inst/testdata/input/testdata.rds")

# Convert to SCE
test_sce <- Seurat::as.SingleCellExperiment(testdata)
colLabels(test_sce) <- Seurat::Idents(testdata)


test_summ <- summarizeAssayByGroup(test_sce, ids = colLabels(test_sce))
test_summ@colData
test_summ@assays@data$mean # gene mean across cell types
test_summ@assays@data$prop.detected # gene prop
# test_summ@assays@data$num.detected # num times detected



# (global) Avg Expr by gene
scuttle::calculateAverage(test_sce)


test_t <- scran::findMarkers(test_sce,
                             groups = colLabels(test_sce),
                             direction = "any",
                             full.stats = TRUE,
                             test.type = "t")
test_t@listData


test_t <- scran::pairwiseTTests(test_sce,
                                groups = colLabels(test_sce))



# Get all possible combinations
pairs <- combn(unique(as.character(colLabels(test_sce))), 2) %>%
    t %>%
    as_tibble()

# t-test by cluster (all genes together)
test_summ@assays@data$mean %>%
    as_tibble() %>%
    select("B", "NK") %>%
    pivot_longer(names_to = "celltype",
                 cols = c("B", "NK"),
                 values_to = "counts")  %>%
    t_test(counts ~ celltype,
           order = c("B", "NK"))


