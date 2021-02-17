# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat.rds")

# Get Omni Resrouces
source("scripts/utils/get_omnipath.R")
omni_resources <- get_omni_resources()

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


# Prepare for ROC
roc_res <- squidpy_sub %>%
    filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(df) calc_curve(df, ground))) # get roc


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
                            seurat_object = call_cellchat,
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



# 3. NATMI
