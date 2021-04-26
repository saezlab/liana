# Load results
spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       # "pval"=FALSE,
                                       "prob"=TRUE
                                       )),
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



# 1. Combine all binary results into heatmap
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


get_BigHeat(top_lists$top_250,
            display_numbers = FALSE,
            silent = FALSE,
            show_rownames = FALSE,
            show_colnames = FALSE,
            legend_breaks = 0:1,
            fontsize = 17,
            drop_levels = TRUE,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            # color = c("gray15", "darkslategray2"),
            border_color = NA,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            treeheight_row = 0,
            treeheight_col = 100)


# 2. Activity by Cell Type Heatmap (Source and Target)
get_activecell(top_lists$top_250)


# Supplementary figures =====
# 3. UpSet Plots  by Tool
names(top_lists$top_250) %>%
    map(function(m_name) top_lists$top_250[[m_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("output/crc_res/plots/upset_tools/{m_name}_upset.png")))


# 4. Upset Plots by Resource
# assign CellPhoneDB to Squidpy default
top250_resource_tool <- get_swapped_list(top_lists$top_250)

# Plot and Save Upsets
names(top250_resource_tool) %>%
    map(function(r_name)
        top250_resource_tool[[r_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("output/crc_res/plots/upset_resources/{r_name}_upset.png"))
    )


# 5. PCA by Rank Frequencies
rank_frequencies <- spec_list %>%
    get_rank_frequencies()
plot_freq_pca(rank_frequencies)


# 6. Get Numbers per Cell Type
get_cellnum("input/crc_data/crc_korean_form.rds")




# Housekeeping measures
housekeep_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       # "pval"=FALSE,
                                       "prob"=TRUE
                                       )),
                  "Connectome" =
                      methods::new("MethodSpecifics",
                                   method_name="Connectome",
                                   method_results = readRDS("output/crc_res/conn_results.rds"),
                                   method_scores=list(
                                       # "weight_sc"=TRUE,
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
                                       # "edge_specificity"=TRUE,
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
                                       # "pvalue"= FALSE,
                                       "means"= TRUE
                                   ))
)

top_housekeep <- get_top_hits(housekeep_list,
                          n_ints=c(250)
)

# 7. Binary Housekeeping heatmap
get_BigHeat(top_housekeep$top_250,
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

# 8. Housekeeping Activity by Cell Type Heatmap
get_activecell(top_housekeep$top_250)


# 9. Housekeeping Rank Avg
housekeep_frequencies <- housekeep_list %>%
    get_rank_frequencies()
plot_freq_pca(housekeep_frequencies)



# Check CellChat P-values
cellchat_alone <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       "pval"=FALSE,
                                       "prob"=TRUE
                                   ))
)

cellchat_hits <- get_top_hits(cellchat_alone,
                          n_ints=c(250)
)

# Check SCA above threshold
sca_hits <- spec_list$SCA@method_results %>%
    map(function(resource) resource %>%
            filter(LRscore >= 0.5))

