# Load results
spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       #"pval"=FALSE,
                                       "prob"=TRUE)),
                  "Squidpy" =
                      methods::new("MethodSpecifics",
                                   method_name="Squidpy",
                                   method_results = readRDS("output/crc_res/squidpy_results.rds"),
                                   method_scores=list(
                                       # "means"=TRUE,
                                       "pvalue"=FALSE
                                       )),
                  "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       "edge_avg_expr"=TRUE,
                                       "edge_specificity"=TRUE)),
                  "iTALK" =
                      methods::new("MethodSpecifics",
                                   method_name="iTALK",
                                   method_results = readRDS("output/crc_res/italk_results.rds"),
                                   method_scores=list(
                                       "weight_comb"=TRUE
                                       )),
                  "Connectome" =
                      methods::new("MethodSpecifics",
                                   method_name="Connectome",
                                   method_results = readRDS("output/crc_res/conn_results.rds"),
                                   method_scores=list(
                                       "weight_sc"=TRUE,
                                       "weight_norm"=TRUE
                                       )),
                  "SCA" = methods::new("MethodSpecifics",
                                       method_name="SCA",
                                       method_results = readRDS("output/crc_res/sca_results.rds"),
                                       method_scores=list(
                                           "LRscore"=TRUE
                                           ))
                  )


# I. Overlap
# Top X Top Hits for each tool
top_lists <- get_top_hits(spec_list, n_ints=c(50, 200, 1000, 5000))


# 1. UpSet Plots and Heatmaps by Tool
names(top_lists$top_1000) %>%
    map(function(m_name) top_lists$top_1000[[m_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("~/Repos/ligrec_decoupleR/output/crc_res/plots/upset_tools/{m_name}_upset.png")))



# 2. Combine all binary results into heatmap
binary_heatm <- get_BigHeat(top_lists$top_1000,
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
                            treeheight_row = 0,
                            treeheight_col = 100)


# 3. Binary PCA
plot_freq_pca(top_lists$top_1000 %>%
                  get_binary_frequencies())


# 4. Upset Plots by Resource
# assign CellPhoneDB to Squidpy default
sig_list_resource <- get_swapped_list(top_lists$top_1000)

# Plot and Save Upsets
names(sig_list_resource) %>%
    map(function(r_name)
        sig_list_resource[[r_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("~/Repos/ligrec_decoupleR/output/crc_res/plots/upset_resources/{r_name}_upset.png"))
    )



# II All Results ranked -------------------------------------------------------
rank_frequencies <- spec_list %>%
    get_rank_frequencies()

# 5. PCA by Rank Frequencies
plot_freq_pca(rank_frequencies)
rank_frequencies
