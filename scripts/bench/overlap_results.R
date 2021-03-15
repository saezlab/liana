# Load Prerequisites
library(pheatmap)
library(jaccard)
library(UpSetR)
library(jaccard)
library(tidyverse)
library(RColorBrewer)
# library(viridis)

# Load results
squidpy_results <- readRDS("output/benchmark/main_run/squidpy_full.rds") %>%
  map(function(ll) ll %>%
        mutate(source = as.character(str_glue("c{source}")),
               target = as.character(str_glue("c{target}")))
      )# issue with clust names
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
        res %>%
        filter(pval <= 0.01) %>%
        mutate(prank = percent_rank(dplyr::desc(prob))) %>%
        filter(prank <= 0.05) %>%
            as_tibble()})

squidpy_sig <- squidpy_results %>%
    map(function(res){
        res %>%
            filter(pvalue <= 0.05) %>%
            as_tibble()})

natmi_sig <- natmi_results %>%
    map(function(res){
        res %>%
        mutate(prank = percent_rank(dplyr::desc(edge_specificity))) %>%
        filter(prank <= 0.005) %>%
            as_tibble()
        })

sca_sig <- sca_results %>%
    map(function(res){
    res %>%
        filter(LRscore >= 0.5) %>% # this is the threshold that they use when they compare
        as_tibble()
        })

italk_sig <- italk_results %>%
    map(function(res){
        res %>%
            mutate(weight_comb = weight_from * weight_to) %>%
            mutate(prank = percent_rank(dplyr::desc(weight_comb))) %>%
            filter(prank < 0.01) %>%
            as_tibble()
    })

conn_sig <- conn_results %>%
    map(function(res){
        res %>% FormatConnectome(max.p = 0.05,
                                 remove.na = TRUE)  %>%
        mutate(prank = percent_rank(desc(weight_sc))) %>%
        filter(prank <= 0.03) %>%
            as_tibble()
    })



# I. Overlap

# 1. UpSet Plots and Heatmaps by Tool
# Significant Hits for each tool
# using specificity measures wheere available (and ranking if necessary)
sig_list <- list("CellChat" = cellchat_sig,
                 "Squidpy" = squidpy_sig,
                 "NATMI" = natmi_sig,
                 "SCA" = sca_sig,
                 "Connectome" = conn_sig,
                 "iTALK" = italk_sig)


# binarize significant results and bind to the same matrix for each tool
# then generate UpsetPlots
names(sig_list) %>%
  map(function(m_name) sig_list[[m_name]] %>%
        prepForUpset() %>%
        plotSaveUset(str_glue("output/benchmark/overlap_plots/upset_tools/{m_name}_upset.png")))




# 2. Combine all binary results into heatmap
binary_heatm <- get_BigHeat(sig_list,
                            display_numbers = FALSE,
                            silent = FALSE,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            legend_breaks = 0:1,
                            fontsize = 18,
                            drop_levels = TRUE,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            color = c("gray15", "darkslategray2"),
                            border_color = NA,
                            clustering_distance_rows = "binary",
                            clustering_distance_cols = "binary",
                            treeheight_row = 0)



# 3. Upset Plots by Resource
# assign CellPhoneDB to Squidpy default
sig_list_resource <- get_swapped_list(sig_list)

# Plot and Save Upsets
names(sig_list_resource) %>%
  map(function(r_name)
    sig_list_resource[[r_name]] %>%
        prepForUpset() %>%
      plotSaveUset(str_glue("output/benchmark/overlap_plots/upset_resources/{r_name}_upset.png"))
  )



# 4. Sig/Top Hits per Cell Pairs PCA
plot_freq_pca(sig_list)


# 5. Combine all (non-binary)
comb <- list("CellChat" = cellchat_results,
             "Squidpy" = squidpy_results,
             "NATMI" = natmi_results,
             "iTALK" = italk_results,
             "SCA" = sca_results,
             "Connectome" = conn_results %>%
               map(function(res){
                 res %>%
                   FormatConnectome(., remove.na = TRUE) %>%
                   as_tibble()
               }))



lnames <- map(names(comb), function(l_name){
  map(names(comb[[l_name]]), function(r_name){
    str_glue("{l_name}_{r_name}")
  })
}) %>% unlist()

comb <- comb %>% purrr::flatten() %>%
  setNames(lnames)
names(comb)


tmp <- map(names(comb),
           function(l_name){
             comb[[l_name]] %>%
               as.data.frame() %>%
               mutate(!!l_name := (.[, 5] - mean(.[, 5])) / sd(.[, 5])) %>%
               select(1:4, !!l_name) %>%
               unite("interaction", source, target, ligand, receptor, sep="_") %>%
               distinct()
}) %>% reduce(., full_join) %>%
  mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
  mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, .)) %>%
  as.data.frame()
tmp

head(tmp)

tmp <- tmp %>% column_to_rownames("interaction") %>% select(3:15)
pca(tmp, labels=colnames(tmp),legendtextsize = 10,axistextsize = 10,dotsize=2)
tmp






