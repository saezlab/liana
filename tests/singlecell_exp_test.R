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



<<<<<<< HEAD
# t-test by cluster (all genes together)
=======
>>>>>>> master
test_summ@assays@data$mean %>%
    as_tibble() %>%
    select("B", "NK") %>%
    pivot_longer(names_to = "celltype",
                 cols = c("B", "NK"),
                 values_to = "counts")  %>%
    t_test(counts ~ celltype,
           order = c("B", "NK"))



<<<<<<< HEAD
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


#
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

seurat_object@meta.data <- seurat_object@meta.data %>%
    rownames_to_column(var = "Bar_Code") %>%
    as_tibble() %>%
    group_by(seurat_annotations) %>%
    slice_sample(prop = 0.1) %>%
    ungroup() %>%
    as.data.frame() %>%
    column_to_rownames("Bar_Code") %>%
    mutate(seurat_annotations = as.factor(as.numeric(seurat_annotations)))

seurat_object <- subset(seurat_object,
                        cells = rownames(seurat_object@meta.data)) %>%
    Seurat::NormalizeData()
length(Seurat::Idents(seurat_object)) # 9 xD

Seurat::Idents(seurat_object) <- seurat_object@meta.data$seurat_annotations

liana_res <- liana_wrap(seurat_object)

liana_res %>% liana_aggregate()
=======
# test CellChat with Mouse ----
ss_liana <- readRDS("~/Downloads/ss_liana.rds")

#' Basic function to convert human to mouse genesymbols (temporary solution)
#' @param op_resource omnipath_resource as obtained via `liana::select_resource`
#'
#' @details adapted from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convert_to_murine <- function(op_resource){
    require("biomaRt")
    require(liana)
    require(tidyverse)

    # query biomaRt databases
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    # obtain tibble with human and murine genesymbol
    symbols_tibble <- getLDS(attributes = c("hgnc_symbol"),
                             filters = "hgnc_symbol",
                             values = union(op_resource$source_genesymbol,
                                            op_resource$target_genesymbol),
                             mart = human,
                             martL = mouse,
                             attributesL = c("mgi_symbol")) %>%
        dplyr::rename(human_symbol = HGNC.symbol,
                      murine_symbol = MGI.symbol) %>%
        as_tibble()

    # intentionally we introduce duplicates, if needed
    # these should be resolved when LIANA is called
    # as inappropriately matched genes will not be assigned any values
    op_resource %>%
        left_join(symbols_tibble, by=c("target_genesymbol"="human_symbol")) %>%
        mutate(target_genesymbol = murine_symbol, .keep = "unused") %>%
        left_join(symbols_tibble, by=c("source_genesymbol"="human_symbol")) %>%
        mutate(source_genesymbol = murine_symbol, .keep = "unused") %>%
        filter(!is.na(target_genesymbol) | !is.na(source_genesymbol)) %>%
        filter(!is.na(target_genesymbol)) %>%
        filter(!is.na(source_genesymbol))
}


require(liana)
require(tidyverse)
require(magrittr)
seurat_object <- readRDS("~/Downloads/ss_liana.rds") %>%
    Seurat::NormalizeData()
op_resource <- select_resource("OmniPath")[[1]] %>%
    convert_to_murine()

liana_res <- liana_wrap(seurat_object,
                        cellchat.params=list(organism="mouse"),
                 resource = "custom",
                 external_resource = op_resource)

# squidpy converts genes to upper, revert this to title (i.e. murine)
liana_res$squidpy %<>%
    mutate_at(.vars = c("ligand", "receptor"), str_to_title)

liana_res %<>% liana_aggregate


>>>>>>> master
