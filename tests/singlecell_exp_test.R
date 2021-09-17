seurat_object <- readRDS("../ligrec_decouple/data/input/cmbc_seurat_test.RDS")
require(SingleCellExperiment)
require(Seurat)

testdata <- readRDS("data/input/test_data.rds")
testdata %<>%
    FindVariableFeatures() %>%
    NormalizeData() %>%
    ScaleData() %>%
    RunPCA(verbose = TRUE) %>%
    FindNeighbors(reduction = "pca") %>%
    FindClusters(resolution = 0.4, verbose = TRUE)


liana_wrap(testdata,
           method = "squidpy")


# convert to singlecell object
sce <- SingleCellExperiment::SingleCellExperiment(
    assays=list(counts = GetAssayData(seurat_object, assay = "ADT", slot = "data")),
    colData=DataFrame(label=seurat_object@meta.data$seurat_clusters)
)


# get OP and filter by ADTs to reduce comp time
op_resource <- select_resource("OmniPath")[[1]] %>%
    filter(source_genesymbol %in% rownames(sce))

liana_res <- liana_wrap(seurat_object,
                        squidpy.params=list(cluster_key = "seurat_clusters"),
                        expr_prop = 0.1,
                        resource = "custom",
                        external_resource = op_resource
)




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


op_resource <- liana::select_resource("OmniPath")






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


# Test ccr data
testdata <- readRDS("data/input/test_data.rds")

xx <- liana_wrap(testdata,
                 method = c("natmi", "connectome", "logfc",
                            "cellchat", "sca", "squidpy"),
                 assay.type = "counts",
                 resource = "CellPhoneDB",
                 parallelize = TRUE,
                 workers = 4)
xx_aggr <- xx %>% liana_aggregate()

xx_aggr %>%
    filter(squidpy.pvalue <= 0.05)
