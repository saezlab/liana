# Robustness Results Analysis
#' @author Daniel Dimitrov, Saezlab

# load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

source("scripts/utils/bench_robust.R")
source("scripts/utils/robust_utils.R")
# sapply(list.files("scripts/pipes/", pattern = ".R", full.names = TRUE), source)


#
# # Load Omni Resources
# omni_resources <-readRDS("input/omni_resources.rds")
#
# # Define Subsampling
# subsampling <- c(1, 0.8, 0.6, 0.4)


# 1. Load Formatted Subsampling Results and generate ROC

# Squidpy
squidpy_sub <- readRDS("output/benchmark/squidpy_sub.rds")
squidpy_roc <- squidpy_sub %>%
    robust_get_roc(predictor_metric = "pvalue", predictor_thresh = 0.05)

# CellChat
cellchat_sub <- readRDS("output/benchmark/cellchat_sub.rds")
cellchat_roc <- cellchat_sub %>%
    robust_get_roc(predictor_metric = "pval", predictor_thresh = 0.05)

# NATMI
natmi_sub <- readRDS("output/benchmark/natmi_sub.rds")
natmi_roc <- natmi_sub %>%
    robust_get_roc(predictor_metric = "edge_avg_expr",
                   predictor_thresh = 0.01,
                   .rank = TRUE)

# Connectome
conn_sub <- readRDS("output/benchmark/connectome_sub.rds")
conn_roc <- conn_sub %>%
    robust_get_roc(predictor_metric = "weight_norm",
                   predictor_thresh = 0.05,
                   .rank = TRUE)

# iTALK
italk_sub <- readRDS("output/benchmark/italk_sub.rds")
italk_roc <- italk_sub %>%
    robust_get_roc(predictor_metric = "weight_comb",
                   predictor_thresh = 0.01,
                   .rank = TRUE)

# SCA
sca_sub <- readRDS("output/benchmark/sca_sub.rds")
sca_roc <- sca_sub %>%
    robust_get_roc(predictor_metric = "LRscore",
                   predictor_thresh = 0.1,
                   .rank = TRUE)



# 2. ROC plots -------------------------------------------------------------
robust_roc_plot(squidpy_roc)
robust_roc_plot(cellchat_roc)
robust_roc_plot(natmi_roc)
robust_roc_plot(conn_roc)
robust_roc_plot(italk_roc)
robust_roc_plot(sca_roc)



# 3. Combine Results into Heatmap
squidpy_heat <- squidpy_roc %>%
    select(resource, subsample, roc) %>%
    ungroup()  %>%
    mutate(resource = str_replace(resource, "CellPhoneDB", "Default")) %>%
    mutate(subsample = str_replace(subsample, "\\.", ",")) %>%
    mutate(alg = "squidpy") %>%
    unite(col = "key", alg, subsample)
squidpy_heat

natmi_heat <- natmi_roc %>%
    select(resource, subsample, roc) %>%
    ungroup() %>%
    mutate(resource = str_replace(resource, "lrc2p", "Default")) %>%
    mutate(alg = "natmi") %>%
    unite(col = "key", alg, subsample)

cellchat_heat <- cellchat_roc %>%
    select(resource, subsample, roc) %>%
    ungroup() %>%
    mutate(alg = "cellchat") %>%
    unite(col = "key", alg, subsample)

conn_roc

conn_heat <- conn_roc %>%
    mutate(name = str_replace(name, "_", "xx")) %>%
    separate(name, into=c("resource", "subsample"), sep = "xx") %>%
    mutate(subsample = str_replace(subsample, "\\.", ",")) %>%
    mutate(alg = "conn") %>%
    unite(col = "key", alg, subsample) %>%
    select(key, roc, resource)



natmi_heat
squidpy_heat
cellchat_heat
conn_heat

heat_data <- readRDS("output/benchmark/heat_data.rds")
# heat_data <- bind_rows(squidpy_heat, cellchat_heat, natmi_heat)
heat_data <- bind_rows(heat_data, conn_heat)
# saveRDS(heat_data ,"output/benchmark/heat_data.rds")


library(pheatmap)
heatp <- heat_data %>%
    filter(!str_detect(key, "subsamp_1")) %>%
    unnest(roc) %>%
    select(key, auc, resource) %>%
    distinct() %>%
    pivot_wider(names_from = resource, values_from = auc) %>%
    column_to_rownames(var = "key") %>%
    pheatmap(.,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = TRUE,
             silent = TRUE)

heatp


squidpy_roc


