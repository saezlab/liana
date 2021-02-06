#### Explore Pipe results

# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

# library(SeuratData


# load python repl
reticulate::repl_python()

# load test seurat object
seurat_object <- readRDS("input/pbmc3k_processed.rds")
table(Idents(seurat_object))

# Get Omni Resrouces
source("./R/utils/get_omnipath.R")
omni_resources <- get_omni_resources()


# 1. CellChat ---------------------------------------------------------------
# Call cellchat and iterate over omni resources
cellchat_results <- omni_resources %>%
    map(function(x) call_cellchat(x,
                                  seurat_object,
                                  thresh = 0.05)) %>%
    setNames(names(omni_resources))


cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = seurat_object,
                                  use_omni_resource = FALSE,
                                  exclude_anns = c("ECM-Receptor",
                                                   "Cell-Cell Contact"),
                                  thresh = 0.05)

# Conclusions:
# Good idea, takes multiple things into considerations that other packages ignore
# cofactors, etc, imputes via PPI, etc.
# However, very difficult to extend Resource, Slow.
# When using OmniPath CellChatDB, we get a slightly different number of
# "significant" hits than when using the original CellChatDB.
# This is liekly due to difference in number of interactions after filtering by annotation
# I opened a Github issue about this and confirmed that what I've done is correct.


# 2.Connectome
connectome_results <- omni_resources %>%
    map(function(op_resource){
        conn <- call_connectome(op_resource,
                                seurat_object,
                                # optional args passed to createConnectom
                                LR.database = 'custom',
                                min.cells.per.ident = 75,
                                p.values = FALSE,
                                calculate.DOR = FALSE)
    })  %>%
    setNames(names(omni_resources)) %>%
    map(function(conn_res) conn_res %>%
            FilterConnectome(.,
                             min.pct = 0.1,
                             min.z = 0.25,
                             remove.na = TRUE))


connectome_default <- call_connectome_default(seurat_object = seurat_object,
                                              min.cells.per.ident = 75,
                                              p.values = FALSE,
                                              calculate.DOR = FALSE) %>%
    FilterConnectome(.,
                     min.pct = 0.1,
                     min.z = 0.25,
                     remove.na = TRUE)



# 3. NATMI ---------------------------------------------------------------------
# Extract data from Seurat Object
write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object,
                                            assay = "RNA",
                                            slot = "data"))) - 1),
          file = "input/test_em.csv",
          row.names = TRUE)
write.csv(Idents(object = seurat_object)  %>%
              enframe(name="barcode", value="annotation"),
          file = "input/test_metadata.csv",
          row.names = FALSE)

# save OmniPath Resource to NATMI format
omni_to_NATMI(omni_resources,
              omni_path = "input/omnipath_NATMI")


# call NATMI
py_set_seed(1004)
natmi_results <- call_natmi(omni_resources,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/test_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/test_metadata.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/NATMI_test")



# 4. Squidpy -------------------------------------------------------------------
library(reticulate)
library(tidyverse)
reticulate::use_python("/home/dbdimitrov/anaconda3/bin/python")
reticulate::source_python("R/pipes/squidpy_source.py")
py$pd <- reticulate::import("pandas")

# prep data for transfer
exprs <- GetAssayData(seurat_object)
meta <- seurat_object[[]]
feature_meta <- GetAssay(seurat_object)[[]]
embedding <- Embeddings(seurat_object, "umap")

# Get Omni Resrouces
source("./R/utils/get_omnipath.R")

# Call Squidpy
reticulate::source_python("R/pipes/squidpy_source.py")
py_set_seed(1004)
py$squidpy_results <- py$call_squidpy(names(omni_resources)[-1], # -CellChatDB returns None
                                      exprs,
                                      meta,
                                      feature_meta,
                                      embedding)





squidpy_pvalues <- py$squidpy_results$pvalues %>% setNames(names(omni_resources)[-1]) #*
squidpy_means <- py$squidpy_results$means %>% setNames(names(omni_resources)[-1]) #*


# reformat squidpy function r
squidpy_reformat <- function(.name,
                             .pval_list,
                             .mean_list){
    x_pval <- .pval_list[[.name]] %>%
        py_to_r() %>%
        pivot_longer(cols = 3:ncol(.),
                     values_to="pvalue",
                     names_to="pair")

    x_mean <- .mean_list[[.name]] %>%
        py_to_r() %>%
        pivot_longer(cols = 3:ncol(.),
                     values_to = "means",
                     names_to = "pair")

    return(left_join(x_pval, x_mean))
}



squidpy_results <- map(names(omni_resources)[-1], #*
                       function(x)
                           squidpy_reformat(.name=x,
                                            .pval_list = squidpy_pvalues,
                                            .mean_list = squidpy_means)) %>%
    setNames(names(omni_resources)[-1]) #*
