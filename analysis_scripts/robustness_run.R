# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

source("scripts/utils/bench_robust.R")
source("scripts/utils/robust_roc.R")
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

# Define Subsampling
subsampling <- c(1, 0.8, 0.6, 0.4)


# 1. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")

# Omni x  Ramilowski x CellPhoneDB
db_list_squidy <- list("OmniPath" = omni_resources$OmniPath,
                       "Ramilowski2015" = omni_resources$Ramilowski2015,
                       "CellPhoneDB" = omni_resources$CellChatDB)

# convert labels to factor (SquidPy)
# breast_cancer@meta.data$lt_id <- as.factor(breast_cancer@meta.data$lt_id)
squidpy_res <- bench_robust(subsampling,
                            call_squidpyR,
                            seurat_object = breast_cancer,
                            omni_resources= db_list_squidy,
                            python_path = "/home/dbdimitrov/anaconda3/bin/python",
                            ident = "seurat_clusters")
# saveRDS(squidpy_res, "output/benchmark/squidpy_res.rds")
squidpy_sub <- squidpy_res %>%
    robust_format_res(., db_list_squidy, subsampling, .res_order = FALSE)
# saveRDS(squidpy_sub, "output/benchmark/squidpy_sub.rds")


# 2. CellChat ------------------------------------------------------------------
source("scripts/pipes/cellchat_pipe.R")

# DBs to Bench
db_list_cellchat <- list("OmniPath" = omni_resources$OmniPath,
                         "Ramilowski2015" = omni_resources$Ramilowski2015,
                         "Default" = NULL)

cellchat_res <- db_list_cellchat %>%
    map(function(db){
        bench_robust(subsampling,
                     lr_call = call_cellchat,
                     op_resource = db,
                     seurat_object = breast_cancer,
                     exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                     nboot = 100,
                     thresh = 1,
                     assay = "SCT")
    })
# saveRDS(cellchat_res, "output/benchmark/cellchat_res.rds")
cellchat_res <- readRDS("output/benchmark/cellchat_res.rds")

cellchat_sub <- cellchat_res %>%
    robust_format_res(., db_list_cellchat, subsampling)
# saveRDS(cellchat_sub, "output/benchmark/cellchat_sub.rds")


# 3. NATMI ---------------------------------------------------------------------
# call NATMI
source("scripts/pipes/NATMI_pipe.R")
reticulate::repl_python()
py_set_seed(1004)

# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources,
#               omni_path = "input/omnipath_NATMI")

db_list_natmi <- list("OmniPath" = omni_resources$OmniPath,
                      "Ramilowski2015" = omni_resources$Ramilowski2015,
                      "Default" = NULL)

natmi_res <- bench_robust(subsampling,
                          lr_call = call_natmi,
                          omni_resources= db_list_natmi,
                          seurat_object = breast_cancer,
                          omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                          natmi_path = "~/Repos/NATMI",
                          em_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_em.csv",
                          ann_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_ann.csv",
                          output_path = "~/Repos/ligrec_decoupleR/output/benchmark/natmi",
                          .write_data = TRUE,
                          .subsampling_pipe = TRUE)
# saveRDS(natmi_res, "output/benchmark/natmi_res.rds")

# load and remove not literature annotated default DB
natmi_res <- readRDS("output/benchmark/natmi_res.rds") %>%
    map(function(x) x %>% purrr::list_modify("lrc2a" = NULL))

# get NATMI res
natmi_sub <- natmi_res %>%
    robust_format_res(., list("Default" = NULL) %>% append(db_list_natmi),
                      subsampling, .res_order = FALSE)
# saveRDS(natmi_sub, "output/benchmark/natmi_sub.rds")



# 4. Connectome ----------------------------------------------------------------
# Freezes with exec/do.call

# Connectome fix - does not work with Seurat4 object
# Requires RNA instead of Spatial as raw counts
db_list_conn <- list("Default" = NULL,
                     "OmniPath" = omni_resources$OmniPath,
                     "Ramilowski2015" = omni_resources$Ramilowski2015)

connectome_res <- map(subsampling, function(ss){
        db_list_conn %>% map(function(db)
            call_connectome(seurat_object = seurat_subsample(breast_cancer, subsampling = ss),
                            op_resource = db,
                            min.cells.per.ident = 10,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            .format = FALSE,
                            assay = 'SCT')
        )
})

# saveRDS(connectome_res, "output/benchmark/connectome_res.rds")
connectome_res <- readRDS("output/benchmark/connectome_res.rds")

conn_sub <- connectome_res %>%
    robust_format_res(., db_list_conn,
                      subsampling,
                      .res_order = FALSE) %>%
    # keep only interactions that are potentially relevant
    mutate(lr_res = lr_res %>% map(function(res){
        FormatConnectome(res,
                         min.pct = 0.1,
                         max.p = 0.05)
    }))
# saveRDS(conn_sub, "output/benchmark/connectome_sub.rds")



# 5. iTALK --------------------------------------------------------------------
db_list_italk <- list("OmniPath" = omni_resources$OmniPath,
                      "Ramilowski2015" = omni_resources$Ramilowski2015,
                      "Default" = NULL)


italk_res <- db_list_italk %>%
    map(function(db)
        bench_robust(subsampling = subsampling,
                     lr_call = call_italk,
                     op_resource = db,
                     breast_cancer,
                     assay = 'SCT',
                     .format = FALSE
        ))
# saveRDS(italk_res, "output/benchmark/italk_res.rds")
italk_res <- readRDS("output/benchmark/italk_res.rds")

italk_sub <- italk_res %>%
    robust_format_res(., db_list_italk,
                      subsampling,
                      .res_order = TRUE) %>%
    # Format and Combine Source and Target means
    mutate(lr_res = lr_res %>% map(function(res){
        FormatiTALK(res) %>%
            mutate(weight_comb = weight_to * weight_from)
    }))
# saveRDS(italk_sub, "output/benchmark/italk_sub.rds")






# 6. SCA ----------------------------------------------------------------------
load("input/LRdb.rda")

db_list_sca <- list("OmniPath" = omni_resources$OmniPath,
                    "Ramilowski2015" = omni_resources$Ramilowski2015,
                    "Default" = NULL)

sca_res <- db_list_sca %>%
    map(function(db)
        bench_robust(subsampling = subsampling,
                     lr_call = call_sca,
                     op_resource = db,
                     breast_cancer,
                     assay = 'SCT',
                     .format = TRUE,
                     s.score = 0,
                     logFC = 0.25
        ))
# saveRDS(sca_res, "output/benchmark/sca_res.rds")
sca_res <- readRDS("output/benchmark/sca_res.rds")

sca_sub <- sca_res %>%
    robust_format_res(., db_list_sca,
                      subsampling,
                      .res_order = TRUE)
# saveRDS(sca_sub, "output/benchmark/sca_sub.rds")


