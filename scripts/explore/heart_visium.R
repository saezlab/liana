# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

# load fibrosis slide
fibrosis_seurat <- readRDS("input/processed_visium/157772.rds")

# Get Omni
source("scripts/utils/get_omnipath.R")
omni_resources <- get_omni_resources()

# Take a look
fibrosis_seurat@assays
fibrosis_seurat@meta.data
fibrosis_seurat@images$slice1@coordinates

# Check Assays
Seurat::GetAssay(fibrosis_seurat)
Seurat::DefaultAssay(fibrosis_seurat) <- "SCT"
Idents(fibrosis_seurat) <- fibrosis_seurat@meta.data$lt_id
# convert labels to factor (SquidPy)
fibrosis_seurat@meta.data$lt_id <- as.factor(fibrosis_seurat@meta.data$lt_id)
fibrosis_seurat <-RenameAssays(fibrosis_seurat, "Spatial" = "RNA")




# 1. CellChat ------------------------------------------------------------------
source("scripts/pipes/cellchat_pipe.R")
cellchat_results <- omni_resources %>%
    map(function(x) call_cellchat(x,
                                  fibrosis_seurat,
                                  nboot = 100,
                                  exclude_anns = c(),
                                  thresh = 1)) %>%
    setNames(names(omni_resources))


cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = fibrosis_seurat,
                                  exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                                  nboot = 100,
                                  thresh = 1,
                                  assay = "SCT")

cellchat_results <- append(cellchat_results,
                           list(cellchat_def = cellchat_default))


# 2. NATMI ---------------------------------------------------------------------
# Extract data from Seurat Object
# write.csv(100 * (exp(as.matrix(GetAssayData(object = fibrosis_seurat,
#                                             assay = "SCT",
#                                             slot = "data"))) - 1),
#           file = "input/fibrosis_em.csv",
#           row.names = TRUE)
# write.csv(Idents(object = fibrosis_seurat)  %>%
#               enframe(name="barcode", value="annotation"),
#           file = "input/fibrosis_metadata.csv",
#           row.names = FALSE)


# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources,
#               omni_path = "input/omnipath_NATMI")


# call NATMI
source("scripts/pipes/NATMI_pipe.R")
reticulate::repl_python()
py_set_seed(1004)
natmi_results <- call_natmi(omni_resources,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/fibrosis_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/fibrosis_metadata.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/fibrosis_natmi")


# 3. Connectome ----------------------------------------------------------------
source("scripts/pipes/connectome_pipe.R")
connectome_results <- omni_resources %>%
    map(function(op_resource){
        conn <- call_connectome(op_resource,
                                fibrosis_seurat,
                                # optional args passed to createConnectome
                                LR.database = 'custom',
                                min.cells.per.ident = 1,
                                p.values = TRUE,
                                calculate.DOR = FALSE,
                                assay = 'SCT')
    })  %>%
    setNames(names(omni_resources))

connectome_default <- call_connectome(op_resource = NULL,
                                      seurat_object = fibrosis_seurat,
                                      min.cells.per.ident = 1,
                                      p.values = TRUE,
                                      calculate.DOR = FALSE,
                                      .format = TRUE,
                                      assay = 'SCT')


# 4. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")
# call squidpy
squidpy_results <- call_squidpyR(seurat_object = fibrosis_seurat,
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 ident = "lt_id")
