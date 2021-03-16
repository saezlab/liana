# Robustness Results Analysis
#' Daniel Dimitrov, Saezlab

# load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

source("R/bench_robust.R")
source("R/robust_roc.R")
# sapply(list.files("scripts/pipes/", pattern = ".R", full.names = TRUE), source)


# 1. Load Formatted Subsampling Results and generate ROC

# Squidpy
squidpy_sub <- readRDS("output/benchmark/squidpy_sub.rds")
squidpy_roc <- squidpy_sub %>%
    robust_get_roc(predictor_metric = "pvalue", predictor_thresh = 0.01)

# CellChat
cellchat_sub <- readRDS("output/benchmark/cellchat_sub.rds")
cellchat_roc <- cellchat_sub %>%
    robust_get_roc(predictor_metric = "pval", predictor_thresh = 0.01)

# NATMI
natmi_sub <- readRDS("output/benchmark/natmi_sub.rds")
natmi_roc <- natmi_sub %>%
    robust_get_roc(predictor_metric = "edge_specificity",
                   predictor_thresh = 0.01,
                   .rank = TRUE)

# Connectome
conn_sub <- readRDS("output/benchmark/connectome_sub.rds")
conn_roc <- conn_sub %>%
    robust_get_roc(predictor_metric = "weight_sc",
                   predictor_thresh = 0.1,
                   .rank = TRUE)

# iTALK
italk_sub <- readRDS("output/benchmark/italk_sub.rds")
italk_roc <- italk_sub %>%
    robust_get_roc(predictor_metric = "weight_comb",
                   predictor_thresh = 0.05,
                   .rank = TRUE)

# SCA
sca_sub <- readRDS("output/benchmark/sca_sub.rds")
sca_roc <- sca_sub %>%
    robust_get_roc(predictor_metric = "LRscore",
                   predictor_thresh = 0.1,
                   .rank = TRUE)


# bind to tibble
roc_tib <- tribble(~alg, ~results,
        "Squidpy", squidpy_roc,
        "CellChat", cellchat_roc,
        "NATMI", natmi_roc,
        "Connectome", conn_roc,
        "iTALK", italk_roc,
        "SCSingalR", sca_roc)


# 2. ROC plots -------------------------------------------------------------
roc_tib %>%
    pmap(function(alg, results){
        message(alg)
        file_name = paste(str_glue("output/benchmark/robust_plots/{alg}_roc.png"))
        png(file_name, width = 900, height = 600)
        print(robust_roc_plot(results))
        dev.off()
    })


# 3. Combine Results into Heatmap
library(pheatmap)

roc_heat.d <- roc_tib %>%
    unnest(results) %>%
    mutate(resource = str_replace(resource, "CellPhoneDB", "Default"),
           resource = str_replace(resource, "lrc2p", "Default")) %>%
    unite(col = "key", alg, subsample) %>%
    filter(!str_detect(key, "subsamp_1")) %>%
    unnest(roc) %>%
    select(key, auc, resource) %>%
    distinct() %>%
    pivot_wider(names_from = resource, values_from = auc) %>%
    column_to_rownames(var = "key")

roc_heat.p <- pheatmap(roc_heat.d,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       treeheight_col = 0,
                       treeheight_row = 0,
                       display_numbers = TRUE,
                       silent = TRUE)

file_name = paste(str_glue("output/benchmark/robust_plots/robustness_heat.png"))
png(file_name, width = 1200, height = 900)
print(roc_heat.p)
dev.off()

# Cut heat
roc_heat.cut <- roc_heat.d %>%
    rownames_to_column("alg_sub") %>%
    filter(!str_detect(alg_sub, "iTALK|NATMI|Connectome")) %>%
    column_to_rownames("alg_sub") %>%
    pheatmap(.,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = TRUE,
             fontsize = 24,
             silent = TRUE)
roc_heat.cut

file_name = paste(str_glue("output/benchmark/robust_plots/robustness_cut_heat.png"))
png(file_name, width = 1200, height = 900)
print(roc_heat.cut)
dev.off()
