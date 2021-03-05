# Load Prerequisites
library(tidyverse)
library(UpSetR)
library(pheatmap)
library(jaccard)
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
        filter(edge_specificity > 0.03) %>%
        mutate(prank = percent_rank(dplyr::desc(edge_avg_expr))) %>%
        filter(prank <= 0.1) %>%
            as_tibble()
        })


sca_sig <- sca_results %>%
    map(function(res){
    res %>%
        filter(LRscore >= 0.5) %>%
        as_tibble()
        })


italk_sig <- italk_results %>%
    map(function(res){
        res %>%
            mutate(weight_comb = weight_from * weight_to) %>%
            mutate(prank = percent_rank(dplyr::desc(weight_comb))) %>%
            filter(prank <= 0.01) %>%
            as_tibble()
    })


source("scripts/pipes/connectome_pipe.R")
conn_sig <- conn_results %>%
    map(function(res){
        res %>% FormatConnectome(max.p = 0.05,
                                 remove.na = TRUE,
                                 min.z = 0.5) %>%
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

# I. Simple Overlap
# 1. UpSet Plots and Heatmaps by Tool
# SquidPy
plotSaveUset(binary_prep$Squidpy,
             "output/benchmark/overlap_plots/squidpy_upset.png")

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

# CellChat
plotSaveUset(binary_prep$CellChat,
             "output/benchmark/overlap_plots/cellchat_upset.png")


# Meaningless - for the BS/Mean score tools


# 2. Upset Plots Each with Default Resource
default_sig <- sig_list %>%
  map(function(tool)
    tool %>% pluck("Default")) %>%
  prepForUpset()

plotSaveUset(default_sig,
             "output/benchmark/overlap_plots/default_sig.png")



# 3. Upset Plots Each tool with OmniPath
omni_sig <- sig_list %>%
  map(function(tool)
    tool %>% pluck("OmniPath")) %>%
  prepForUpset()

plotSaveUset(omni_sig,
             "output/benchmark/overlap_plots/omni_sig.png")

# remove CellChat
sig_list_excl <- sig_list
sig_list_excl$CellChat <- NULL
minus_cell_chat <- sig_list_excl %>%
  map(function(tool)
    tool %>% pluck("OmniPath")) %>%
  prepForUpset()

plotSaveUset(minus_cell_chat,
             "output/benchmark/overlap_plots/omni_minus_cc.png")


# 4. Upset Plots Each tool with Ramilowski
ramilowski_sig <- sig_list %>%
  map(function(tool)
    tool %>% pluck("Ramilowski2015")) %>%
  prepForUpset()

plotSaveUset(ramilowski_sig,
             "output/benchmark/overlap_plots/ramilowski_sig.png")


# 5. Upset Plots Each tool with Random
random_sig <- sig_list %>%
  map(function(tool)
    tool %>%
      pluck("Random")) %>%
  purrr::list_modify("Squidpy" = NULL) %>%
  prepForUpset()


plotSaveUset(random_sig,
             "output/benchmark/overlap_plots/random_sig.png")



# 6. Upset Plots Each tool with CellChatDB
cellchatdb_sig <- sig_list %>%
  map(function(tool)
    tool %>%
      pluck("CellChatDB")) %>%
  # purrr::list_modify("CellChat" = NULL) %>%
  prepForUpset()


plotSaveUset(cellchatdb_sig,
             "output/benchmark/overlap_plots/cellchatdb_sig.png")


# II. Overlap Scaled by Sig. Hits from SquidPy
tmp <- sig_list %>%
  map(function(db)
    db %>%
      pluck("Default")) %>%
  map(function(tool) tool %>%
        select(source, target)) %>%
  enframe() %>%
  mutate(proportions = value %>%
           map(function(cols) cols %>%
                 unite(source, target, col = "pairs"))) %>%
  select(name, proportions)
