# Load Prerequisites
library(tidyverse)
library(UpSetR)
library(pheatmap)
source("scripts/utils/plot_utils.R")


# Load results
squidpy_results <- readRDS("output/benchmark/main_run/squidpy_full.rds") %>%
    plyr::rename(., c("CellPhoneDB" = "Default"))
cellchat_results <- readRDS("output/benchmark/main_run/cellchat_full.rds")
natmi_results <- readRDS("output/benchmark/main_run/natmi_full.rds")  %>%
    purrr::list_modify("lrc2a" = NULL) %>%
    plyr::rename(., c("lrc2p" = "Default"))
sca_results <- readRDS("output/benchmark/main_run/sca_full.rds")
italk_results <- readRDS("output/benchmark/main_run/italk_full.rds")
conn_results <- readRDS("output/benchmark/main_run/conn_full.rds")


# Filter for 'significant' hits and format to Upset
cellchat_sig <- cellchat_results  %>%
    map(function(res){
        res %>% filter(pval < 0.01) %>%
            as_tibble()})


squidpy_sig <- squidpy_results %>%
    map(function(res){
        res %>% mutate(source = str_glue("c{source}")) %>%
            filter(pvalue < 0.05) %>%
            as_tibble()})


natmi_sig <- natmi_results %>%
    map(function(res){
        res %>%
            mutate(prank = percent_rank(edge_avg_expr)) %>%
            filter(prank <= 0.1) %>%
            filter(edge_specificity > 0.05) %>%
            as_tibble()
        })


sca_sig <- sca_results %>%
    map(function(res){
    res %>%
        filter(LRscore >= 0.6) %>%
        as_tibble()
        })


italk_sig <- italk_results %>%
    map(function(res){
        res %>%
            mutate(weight_comb = weight_from * weight_to) %>%
            mutate(prank = percent_rank(weight_comb)) %>%
            filter(prank <= 0.01) %>%
            as_tibble()
    })


source("scripts/pipes/connectome_pipe.R")
conn_sig <- conn_results %>%
    map(function(res){
        res %>% FormatConnectome(max.p = 0.05,
                                 remove.na = TRUE,
                                 min.z = 0.3) %>%
            as_tibble()
    })


# Significant Hits
sig_list <- list("CellChat" = cellchat_sig,
                    "Squidpy" = squidpy_sig,
                    "NATMI" = natmi_sig,
                    "SCA" = sca_sig,
                    "Connectome" = conn_sig,
                    "iTALK" = italk_sig)

binary_prep <- sig_list %>%
  map(function(df) prepForUpset(df))



# 1. UpSet Plots and Heatmaps by Tool
upset(binary_prep$Squidpy, nsets = ncol(binary_prep$Squidpy), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")


heatm_test <- binary_prep$Squidpy %>%
  as_tibble() %>%
  column_to_rownames("interaction") %>%
  pheatmap(.,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           treeheight_col = 0,
           treeheight_row = 0,
           display_numbers = FALSE,
           silent = TRUE,
           show_rownames = FALSE
           )
heatm_test

# Meaningless - for the BS/Mean score tools


# 2. Upset Plots Each with Default Resource
xd <- sig_list %>%
  map(function(tool)
    tool %>% pluck("Default")) %>%
  prepForUpset()


upset(xd, nsets = ncol(binary_prep$Squidpy), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")

tmp <- xd[[-"CellChat"]] %>%
  as_tibble() %>%
  column_to_rownames("interaction") %>%
  pheatmap(.,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           treeheight_col = 0,
           treeheight_row = 0,
           display_numbers = FALSE,
           silent = TRUE,
           show_rownames = FALSE
  )
tmp




# 3. Upset Plots Each tool with OmniPath
xd <- sig_list %>%
  map(function(tool)
    tool %>% pluck("OmniPath")) %>%
  prepForUpset()


upset(xd, nsets = ncol(binary_prep$Squidpy), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")

tmp <- xd[[-"CellChat"]] %>%
  as_tibble() %>%
  column_to_rownames("interaction") %>%
  pheatmap(.,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           treeheight_col = 0,
           treeheight_row = 0,
           display_numbers = FALSE,
           silent = TRUE,
           show_rownames = FALSE
  )
tmp


# 4. Upset Plots Each tool with Random
xd <- sig_list %>%
  map(function(tool)
    tool %>% pluck("Ramilowski2015")) %>%
  prepForUpset()


upset(xd, nsets = ncol(binary_prep$Squidpy), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")

tmp <- xd[[-"CellChat"]] %>%
  as_tibble() %>%
  column_to_rownames("interaction") %>%
  pheatmap(.,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           treeheight_col = 0,
           treeheight_row = 0,
           display_numbers = FALSE,
           silent = TRUE,
           show_rownames = FALSE
  )
tmp


# 5. Heatmap of Each Tool, Each Resource (Sig Hits)


# 6. Heatmap of Each Tool,
