# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

sapply(list.files("scripts/pipes/", pattern = ".R", full.names = TRUE), source)

# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat323.rds")
# Fix for NATMI
clust.anns <- c("c0", "c1","c2",
                "c3","c4","c5",
                "c6", "c7", "c8",
                "c9", "c10", "c11", "c12")
names(clust.anns) <- levels(breast_cancer)
breast_cancer <- RenameIdents(breast_cancer, clust.anns)

# Get Omni Resrouces
# source("scripts/utils/get_omnipath.R")
# omni_resources <- get_omni_resources()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <-readRDS("input/omni_resources.rds")



# Get Random DB
source("scripts/utils/shuffle_omnipath.R")
op_random <- shuffle_omnipath(omni_resources$OmniPath)

# Append shuffled omni db and default to omni_resources
omni_resources_plus <- append(list("Random" = op_random,
                                   "Default" = NULL),
                              omni_resources)


# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = breast_cancer,
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 ident = "seurat_clusters")
# saveRDS(squidpy_results, "output/benchmark/main_run/squidpy_full.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources = omni_resources,
              # omni_path = "input/omnipath_NATMI")
py_set_seed(1004)
natmi_results <- call_natmi(omni_resources = omni_resources_plus,
                            seurat_object = breast_cancer,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_em_subsample_1.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_ann_subsample_1.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/benchmark/natmi_full",
                            .write_data = FALSE,
                            .subsampling_pipe = FALSE
                            )
# saveRDS(natmi_results, "output/benchmark/main_run/natmi_full.rds")


# 3. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources_plus %>%
    map(function(db) call_cellchat(db,
                                  seurat_object = breast_cancer,
                                  nboot = 100,
                                  exclude_anns = c(),
                                  thresh = 1,
                                  assay = "SCT")) %>%
    setNames(names(omni_resources_plus))
# saveRDS(cellchat_results, "output/benchmark/main_run/cellchat_full.rds")


# 4. SCA ----------------------------------------------------------------------
sca_results <- omni_resources_plus %>%
    map(function(db)
        call_sca(op_resource = db,
                 breast_cancer,
                 assay = 'SCT',
                 .format = TRUE,
                 s.score = 0,
                 logFC = 0.25
        ))
# saveRDS(sca_results, "output/benchmark/main_run/sca_full.rds")


# 5. Connectome ----------------------------------------------------------------
conn_results <- omni_resources_plus %>%
    map(function(db)
    call_connectome(seurat_object = breast_cancer,
                    op_resource = db,
                    min.cells.per.ident = 10,
                    p.values = TRUE,
                    calculate.DOR = FALSE,
                    .format = FALSE,
                    assay = 'SCT'))
# saveRDS(conn_results, "output/benchmark/main_run/conn_full.rds")


# 6. iTALK
italk_results <- omni_resources_plus %>%
    map(function(db)
        call_italk(op_resource = db,
                   breast_cancer,
                   assay = 'SCT',
                   .format = TRUE
        ))
# saveRDS(italk_results, "output/benchmark/main_run/italk_full.rds")
