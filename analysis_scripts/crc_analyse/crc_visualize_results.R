# Load results
spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       #"pval"=FALSE,
                                       "prob"=TRUE)),
                  "Connectome" =
                      methods::new("MethodSpecifics",
                                   method_name="Connectome",
                                   method_results = readRDS("output/crc_res/conn_results.rds"),
                                   method_scores=list(
                                       "weight_sc"=TRUE #,
                                       # "weight_norm"=TRUE
                                   )),
                  "iTALK" =
                      methods::new("MethodSpecifics",
                                   method_name="iTALK",
                                   method_results = readRDS("output/crc_res/italk_results.rds"),
                                   method_scores=list(
                                       "weight_comb"=TRUE
                                   )),
                  "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       "edge_specificity"=TRUE # ,
                                       # "edge_avg_expr"=TRUE
                                       )),
                  "SCA" = methods::new("MethodSpecifics",
                                       method_name="SCA",
                                       method_results = readRDS("output/crc_res/sca_results.rds"),
                                       method_scores=list(
                                           "LRscore"=TRUE
                                           )),
                  "Squidpy" =
                      methods::new("MethodSpecifics",
                                   method_name="Squidpy",
                                   method_results = readRDS("output/crc_res/squidpy_results.rds"),
                                   method_scores=list(
                                       # "means"=TRUE,
                                       "pvalue"=FALSE
                                   ))
                  )

# I. Overlap
# Top X Top Hits for each tool
top_lists <- get_top_hits(spec_list,
                          n_ints=c(250)
                          )


# 1. UpSet Plots and Heatmaps by Tool
names(top_lists$top_250) %>%
    map(function(m_name) top_lists$top_250[[m_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("~/Repos/ligrec_decoupleR/output/crc_res/plots/upset_tools/{m_name}_upset.png")))


# 2. Upset Plots by Resource
# assign CellPhoneDB to Squidpy default
top250_resource_tool <- get_swapped_list(top_lists$top_250)

# Plot and Save Upsets
names(top250_resource_tool) %>%
    map(function(r_name)
        top250_resource_tool[[r_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("~/Repos/ligrec_decoupleR/output/crc_res/plots/upset_resources/{r_name}_upset.png"))
    )


# 3. Combine all binary results into heatmap
binary_heatm <- get_BigHeat(top_lists$top_250,
                            display_numbers = FALSE,
                            silent = FALSE,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            legend_breaks = 0:1,
                            fontsize = 17,
                            drop_levels = TRUE,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            color = c("gray15", "darkslategray2"),
                            border_color = NA,
                            clustering_distance_rows = "binary",
                            clustering_distance_cols = "binary",
                            treeheight_row = 0,
                            treeheight_col = 100)

# 4. Activity by Cell Type Heatmap (Source and Target)
get_activecell(top_lists$top_250)



# 5. PCA by Rank Frequencies
rank_frequencies <- spec_list %>%
    get_rank_frequencies()
plot_freq_pca(rank_frequencies)


# 6. Get Numbers per Cell Type
get_cellnum("input/crc_data/crc_korean_form.rds")




# Supp
# Load results
spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       # "pval"=FALSE,
                                       "prob"=TRUE)),
                  "Connectome" =
                      methods::new("MethodSpecifics",
                                   method_name="Connectome",
                                   method_results = readRDS("output/crc_res/conn_results.rds"),
                                   method_scores=list(
                                       "weight_sc"=TRUE,
                                       "weight_norm"=TRUE
                                   )),
                  "iTALK" =
                      methods::new("MethodSpecifics",
                                   method_name="iTALK",
                                   method_results = readRDS("output/crc_res/italk_results.rds"),
                                   method_scores=list(
                                       "weight_comb"=TRUE
                                   )),
                  "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       "edge_specificity"=TRUE,
                                       "edge_avg_expr"=TRUE
                                   )),
                  "SCA" = methods::new("MethodSpecifics",
                                       method_name="SCA",
                                       method_results = readRDS("output/crc_res/sca_results.rds"),
                                       method_scores=list(
                                           "LRscore"=TRUE
                                       )),
                  "Squidpy" =
                      methods::new("MethodSpecifics",
                                   method_name="Squidpy",
                                   method_results = readRDS("output/crc_res/squidpy_results.rds"),
                                   method_scores=list(
                                       "means"=TRUE,
                                       "pvalue"=FALSE
                                   ))
)

top_lists <- get_top_hits(spec_list,
                          n_ints=c(250)
)


# Supp. Check CellChat pva-lues
cellchat_alone <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       "pval"=FALSE,
                                       "prob"=TRUE))
                      )

top_cellchat <- get_top_hits(cellchat_alone,
                             n_ints=c(250)
                             )
