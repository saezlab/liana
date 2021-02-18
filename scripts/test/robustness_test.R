# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

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
# source("scripts/utils/shuffle_omnipath.R")
# op_random <- shuffle_omnipath(omni_resources$OmniPath)

source("scripts/utils/bench_robust.R")
# DBs to Bench
db_list <- list("CellChatDB" = omni_resources$CellChatDB,
                "CellTalkDB" = omni_resources$CellTalkDB)

# Define Subsampling
# subsampling <- c(1, 0.8, 0.6, 0.4)
subsampling <- c(1, 0.75, 0.5)


# 1. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")

# convert labels to factor (SquidPy)
# breast_cancer@meta.data$lt_id <- as.factor(breast_cancer@meta.data$lt_id)
squidpy_res <- bench_robust(subsampling,
             call_squidpyR,
             seurat_object = breast_cancer,
             omni_resources= db_list,
             python_path = "/home/dbdimitrov/anaconda3/bin/python",
             ident = "seurat_clusters")  %>%
    purrr::flatten() %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("{name}_subsamp_{rep(subsampling, each = length(db_list))}"))

squidpy_sub <- squidpy_res %>%
    filter(str_detect(name, "CellChatDB"))
squidpy_sub

# Define Truth using full data set
ground <- squidpy_sub[1,]$lr_res[[1]] %>%
    mutate(truth = if_else(pvalue > 0.05 | is.na(pvalue), 0, 1)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth)

# negative vastly more than positve - need to downsample for PR
summary(as.factor(ground$truth))


# Get ROC
roc_res <- squidpy_sub %>%
    filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(df)
        calc_curve(df, ground, predictor_metric = "pvalue"))) # get roc


ggplot(roc_res %>%
           unnest(roc), aes(x = 1-specificity,
               y = sensitivity,
               colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")



# 2. CellChat ------------------------------------------------------------------
source("scripts/pipes/cellchat_pipe.R")
cellchat_default <- call_cellchat(op_resource = NULL,
                                  seurat_object = breast_cancer,
                                  exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                                  nboot = 100,
                                  thresh = 1,
                                  assay = "SCT")

# DBs to Bench
db_list <- list("Ramilowski2015" = omni_resources$Ramilowski2015)


cellchat_res <- bench_robust(subsampling,
                            lr_call = call_cellchat,
                            op_resource = omni_resources$Ramilowski2015,
                            seurat_object = breast_cancer,
                            exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                            nboot = 100,
                            thresh = 1,
                            assay = "SCT") %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("subsamp_{rep(subsampling, each = length(db_list))}"))


# Define Truth using full data set
ground <- cellchat_res[1,]$lr_res[[1]] %>%
    mutate(truth = if_else(pval > 0.05 | is.na(pval), 0, 1)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth)

# negative vastly more than positve - need to downsample for PR
summary(as.factor(ground$truth))


# Prepare for ROC
roc_res <- cellchat_res %>%
    filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(df)
        calc_curve(df, ground, predictor_metric = "pval"))) # get roc



ggplot(roc_res %>%
           unnest(roc), aes(x = 1-specificity,
                            y = sensitivity,
                            colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")



# 3. NATMI ---------------------------------------------------------------------
# call NATMI
source("scripts/pipes/NATMI_pipe.R")
reticulate::repl_python()
py_set_seed(1004)

# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources,
#               omni_path = "input/omnipath_NATMI")

db_list <- list("CellChatDB" = omni_resources$CellChatDB)

natmi_res <- bench_robust(subsampling,
                          lr_call = call_natmi,
                          omni_resources= db_list,
                          seurat_object = breast_cancer,
                          omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                          natmi_path = "~/Repos/NATMI",
                          em_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_em.csv",
                          ann_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_ann.csv",
                          output_path = "~/Repos/ligrec_decoupleR/output/bc_natmi_subsample",
                          .write_data = TRUE,
                          .subsampling_pipe = TRUE)

natmi_sub <- natmi_res %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("subsamp_{rep(subsampling, each = length(db_list))}"))

xd <- natmi_res %>%
    map(function(x) pull(x)) %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("subsamp_{rep(subsampling, each = length(db_list))}"))

# Define Truth using full data set
ground <- natmi_sub[1,]$lr_res[[1]][[1]] %>%
    mutate(top_ntile = ntile(edge_avg_expr, 100)) %>%
    mutate(truth = if_else(top_ntile > 20 & edge_specificity > 0.05, 1, 0)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth)

# attempt to keep the same proportion of hits as other tools
summary(as.factor(ground$truth))

# Prepare for ROC
roc_res <- natmi_sub %>%
    dplyr::filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(.df)
        calc_curve(.df, ground, predictor_metric = "edge_avg_expr"))) # get roc


ggplot(roc_res %>%
           unnest(roc), aes(x = 1-specificity,
                            y = sensitivity,
                            colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")



# 4. Connectome ----------------------------------------------------------------
source("scripts/pipes/connectome_pipe.R")

# Freezes with exec/do.call

# Connectome fix - does not work with Seurat4 object
# Requires RNA instead of Spatial as raw counts
db_list <- list("Default" = NULL)

seurat_subsets <- map(subsampling, function(ss){
    seurat_object = seurat_subsample(breast_cancer, subsampling = ss)
    })


connectome_res <- seurat_subsets %>%
    map(function(seurat_sub){
        call_connectome(seurat_object = seurat_sub,
                        op_resource = NULL,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = FALSE,
                        assay = 'SCT')
    }) %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("subsamp_{rep(subsampling, each = length(db_list))}"))



# connectome_results <- omni_resources %>%
#     map(function(op_resource){
#         conn <- call_connectome(op_resource,
#                                 breast_cancer,
#                                 # optional args passed to createConnectome
#                                 min.cells.per.ident = 1,
#                                 p.values = TRUE,
#                                 calculate.DOR = FALSE,
#                                 assay = 'SCT')
#     })  %>%
#     setNames(names(omni_resources))


# Freezes and hangs for me on Seurat::ScaleData when called from exec/do.call
# connectome_res <- bench_robust(subsampling = subsampling,
#                                lr_call = call_connectome,
#                                op_resource = omni_resources$CellChatDB,
#                                seurat_object = breast_cancer,
#                                min.cells.per.ident = 1,
#                                p.values = FALSE,
#                                calculate.DOR = FALSE,
#                                .format = FALSE,
#                                assay = 'SCT')

glimpse(connectome_sub[1,]$lr_res[[1]])

ground <- connectome_sub[1,]$lr_res[[1]] %>%
    mutate(top_ntile = ntile(weight_sc, 1000)) %>%
    mutate(truth = if_else((p_val_adj.lig <= 0.05 & p_val_adj.rec <= 0.05) &
                               top_ntile > 995 &
                               (percent.source > 0.1 & percent.target > 0.1),
                           1,
                           0)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth) %>%
    na.omit()

summary(as.factor(ground$truth))


# attempt to keep the same proportion of hits as other tools
summary(as.factor(ground$truth))
connectome_res$lr_res[[2]]

# Prepare for ROC
roc_res <- connectome_res %>%
    dplyr::filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(.df){
        .df %>%
            na.omit() %>%
            calc_curve(., ground, predictor_metric = "weight_sc") # get roc
    }))


xd <- connectome_res$lr_res[[1]]



ggplot(roc_res %>%
           unnest(roc), aes(x = 1-specificity,
                            y = sensitivity,
                            colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")


