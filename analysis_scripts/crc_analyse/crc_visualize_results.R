# To become an .rmd

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
                          n_ints=c(100,
                                   250,
                                   500,
                                   1000)
                          )



# 1. Combine all binary results into heatmap
binary_heatm <- get_BinaryHeat(top_lists$top_500,
                               display_numbers = FALSE,
                               silent = FALSE,
                               show_rownames = FALSE,
                               show_colnames = FALSE,
                               legend_breaks = 0:1,
                               fontsize = 17,
                               drop_levels = TRUE,
                               cluster_rows = FALSE,
                               cluster_cols = TRUE,
                               color = c("gray15", "darkslategray2"),
                               border_color = NA,
                               clustering_distance_rows = "binary",
                               clustering_distance_cols = "binary",
                               treeheight_row = 0,
                               treeheight_col = 100)
binary_heatm


# 2. Activity by Cell Type Heatmap (Source and Target)
get_activecell(top_lists$top_500, cap_value = 0.2)


# Supplementary figures =====
# 3. UpSet Plots  by Tool
names(top_lists$top_500) %>%
    map(function(m_name) top_lists$top_500[[m_name]] %>%
            prepForUpset() %>%
            plotSaveUset(str_glue("output/crc_res/plots/upset_tools/{m_name}_upset.png")))

# 4. Upset Plots by Resource
# assign CellPhoneDB to Squidpy default
top250_resource_tool <- get_swapped_list(top_lists$top_500)

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


# 7. Similarity Heatmap and Stats
get_simdist_heatmap(top_lists$top_500,
                    sim_dist = "simil",
                    method = "Jaccard",
                    diag = TRUE,
                    upper = TRUE)


jaccard_per_mr <- simdist_resmet(top_lists$top_500,
                                 sim_dist = "simil",
                                 method = "Jaccard")
jac <- list_stats(meth = jaccard_per_mr$meth,
                  reso = jaccard_per_mr$reso)



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
                      "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       # "edge_specificity"=TRUE,
                                       "edge_avg_expr"=TRUE
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
                              n_ints=c(500))

# 8. Binary Housekeeping heatmap
get_BinaryHeat(top_housekeep$top_500,
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


# 9. Housekeeping Activity by Cell Type Heatmap
get_activecell(top_housekeep$top_500)


# 10. Housekeeping Rank Avg
housekeep_list %>%
    get_rank_frequencies() %>%
    plot_freq_pca()

# 11. Housekeeping Bray Curtis Info
get_simdist_heatmap(top_housekeep$top_500,
                    sim_dist = "simil",
                    method = "Jaccard",
                    diag = TRUE,
                    upper = TRUE)

jaccard_house <- simdist_resmet(top_housekeep$top_500,
                                 sim_dist = "simil",
                                 method = "Jaccard")
list_stats(meth = jaccard_house$meth,
           reso = jaccard_house$reso)




# 12. Addititional Checks
# Universe of only those 4 methods
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
                  "NATMI" =
                      methods::new("MethodSpecifics",
                                   method_name="NATMI",
                                   method_results = readRDS("output/crc_res/natmi_results.rds"),
                                   method_scores=list(
                                       "edge_specificity"=TRUE # ,
                                       # "edge_avg_expr"=TRUE
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

top_lists <- get_top_hits(spec_list,
                          n_ints=c(500))


jaccard_per_mr <- simdist_resmet(top_lists$top_500,
                                 sim_dist = "simil",
                                 method = "Jaccard")
list_stats(meth = jaccard_per_mr$meth,
           reso = jaccard_per_mr$reso)



# Check CellChat P-values
cc_hits <- spec_list$CellChat@method_results %>%
    map(function(resource) resource %>%
            filter(pval == 0))
cc_hits

# Check SCA above threshold
sca_hits <- spec_list$SCA@method_results %>%
    map(function(resource) resource %>%
            filter(LRscore >= 0.5))
sca_hits

sq_hits <- spec_list$Squidpy@method_results %>%
    map(function(resource) resource %>%
            filter(pvalue <= 0.05))
sq_hits



# Supp Note 2 Complexes
complex_resources <-  c("Baccin2019",
                        "CellChatDB",
                        "CellPhoneDB",
                        "ICELLNET",
                        "Default")

spec_list <- list("CellChat" =
                      methods::new("MethodSpecifics",
                                   method_name="CellChat",
                                   method_results = readRDS("output/crc_res/cellchat_results.rds"),
                                   method_scores=list(
                                       "pval"=FALSE
                                   )),
                  "Squidpy" =
                      methods::new("MethodSpecifics",
                                   method_name="Squidpy",
                                   method_results = readRDS("output/crc_res/squidpy_results.rds"),
                                   method_scores=list(
                                       "pvalue"=FALSE
                                   ))
                  )

spec_list$CellChat@method_results %<>% keep(names(.) %in% complex_resources)
spec_list$Squidpy@method_results %<>% keep(names(.) %in% complex_resources[-5])

# Get Significant hits for Squidpy and CellChat
sig_list <- get_top_hits(spec_list,
                         n_ints=c(5000))

# Sig hit heatmap
get_BinaryHeat(sig_list$top_5000,
               display_numbers = FALSE,
               silent = FALSE,
               show_rownames = FALSE,
               show_colnames = FALSE,
               legend_breaks = 0:1,
               fontsize = 17,
               drop_levels = TRUE,
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               color = c("gray15", "darkslategray2"),
               border_color = NA,
               clustering_distance_rows = "binary",
               clustering_distance_cols = "binary",
               treeheight_row = 0,
               treeheight_col = 100)



# 2. Activity by Cell Type Heatmap (Source and Target)
get_activecell(sig_list$top_5000)


# 3. Percentages of Complexes
compl_perc <- sig_list$top_5000 %>%
    enframe(name = "method") %>%
    unnest(value) %>%
    mutate(resource = names(value)) %>%
    unnest(value) %>%
    select(method, resource, ligand, receptor) %>%
    group_by(method, resource) %>%
    add_count(name = "total") %>%
    filter((str_detect(receptor, "_") | str_detect(ligand, "_"))) %>%
    add_count(name = "complex") %>%
    mutate(prop = complex/total) %>%
    select(resource, method, prop) %>%
    distinct()

mean(compl_perc$prop)


sig_list$top_5000$CellChat$Default %>%
    filter(str_detect(receptor, "_"))
