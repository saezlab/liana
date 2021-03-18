#devtools::install_local('/home/james/sciebo/LR_Benchmark/code/SingleCellSignalR_v1/SingleCellSignalR/')
require(Seurat)
library(SCAomni) ## Install from sciebo

source('/home/james/sciebo/Cell_Cell_Investigation-master/scripts/utils/get_omnipath.R')

seurat_object <- readRDS('/home/james/sciebo/LR_Benchmark/data/pbmc3k_processed.rds')
sel <- dataset$CellChatDB
sca_res = call_sca(sel,seurat_object,assay='RNA',.format=TRUE)

# Second Test
sca_res <- call_sca(op_resource,
                    breast_cancer,
                    assay='SCT',
                    .format=TRUE,
                    s.score = 0,
                    logFC = log2(0.25))


# Default
load("input/LRdb.rda")
sca_default <- call_sca(op_resource = LRdb,
                        breast_cancer,
                        assay='SCT',
                        .format=TRUE,
                        .default_db = TRUE,
                        s.score = 0,
                        logFC = log2(0.25))
