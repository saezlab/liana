# Load results
squidpy_results <- readRDS("output/benchmark/main_run/squidpy_full.rds")
cellchat_results <- readRDS("output/benchmark/main_run/cellchat_full.rds")
natmi_results <- readRDS("output/benchmark/main_run/natmi_full.rds")
sca_results <- readRDS("output/benchmark/main_run/sca_full.rds")
italk_results <- readRDS("output/benchmark/main_run/italk_full.rds")
conn_results <- readRDS("output/benchmark/main_run/conn_full.rds")


# I. Overlap

# Significant/Top Hits for each tool
squidpy_sig <- squidpy_results %>%
  map(function(res){
    res %>%
      filter(pvalue <= 0.05) %>%
      as_tibble()
  })

cellchat_sig <- cellchat_results %>%
  map(function(res){
    res %>%
      mutate(pval = p.adjust(pval, method = "BH")) %>%
      filter(pval <= 0.00) %>%
      mutate(prank = percent_rank(dplyr::desc(prob))) %>%
      filter(prank <= 0.1) %>%
      as_tibble()})

natmi_sig <- natmi_results %>%
  map(function(res){
    res %>%
      filter(source != target) %>%
      mutate(prank = percent_rank(dplyr::desc(edge_specificity))) %>%
      filter(prank <= 0.01) %>%
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
      filter(qval_from <= 0.05 & qval_to <= 0.05) %>%
      as_tibble()
  })

conn_sig <- conn_results %>%
  map(function(res){
    res %>%
      filter(p_val_adj.lig <= 0.05 & p_val_adj.rec <= 0.05) %>%
      mutate(prank = percent_rank(desc(weight_sc))) %>%
      filter(prank <= 0.1) %>%
      as_tibble()
  })



sig_list <- list("CellChat" = cellchat_sig,
                 "Squidpy" = squidpy_sig,
                 "NATMI" = natmi_sig,
                 "iTALK" = italk_sig,
                 "Connectome" = conn_sig,
                 "SCA" = sca_sig) # order for hm



# 1. UpSet Plots and Heatmaps by Tool
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


# 4. Binary PCA
plot_freq_pca(sig_list %>%
                get_binary_frequencies())

# 4.1 Binary PCA - SCA
plot_freq_pca(sig_list %>%
                purrr::list_modify("SCA" = NULL) %>%
                get_binary_frequencies())



# II All Results ranked -------------------------------------------------------
squidpy_full <- squidpy_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="pvalue",
                           .desc_order = FALSE)
  )

# NATMI
natmi_full <- natmi_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="edge_specificity",
                           .desc_order = TRUE)
  )

# SCA
sca_full <- sca_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="LRscore",
                           .desc_order = TRUE)
  )

# Connectome
conn_full <- conn_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="weight_sc",
                           .desc_order = TRUE)
  )

# CellChat
cellchat_full <- cellchat_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="pval",
                           .desc_order = FALSE)
  )

# iTALK
italk_full <- italk_results %>%
  map(function(res)
    res %>%
      format_rank_frequencies(score_col="weight_comb",
                           .desc_order = TRUE)
  )


# Combine all into list and get frequencies per rank
rank_frequencies <- (list("CellChat" = cellchat_full,
                   "Squidpy" = squidpy_full,
                   "NATMI" = natmi_full,
                   "iTALK" = italk_full,
                   "Connectome" = conn_full,
                   "SCA" = sca_full)) %>%
  get_rank_frequencies()

# 5. PCA by Rank Frequencies
plot_freq_pca(rank_frequencies)
rank_frequencies


