# LM/Corr with NES fig + bar plots
# Combine all into list and get frequencies per rank
# II All Results ranked -------------------------------------------------------
squidpy_p_full <- squidpy_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="pvalue",
                                    .desc_order = FALSE)
    )

squidpy_m_full <- squidpy_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="means",
                                    .desc_order = TRUE)
    )

# NATMI
natmi_sw_full <- natmi_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="edge_specificity",
                                    .desc_order = TRUE)
    )


natmi_epw_full <- natmi_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="edge_avg_expr",
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
conn_sw_full <- conn_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="weight_sc",
                                    .desc_order = TRUE)
    )

conn_epw_full <- conn_results %>%
    map(function(res)
        res %>%
            format_rank_frequencies(score_col="weight_norm",
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
                          "Squidpy.pvals" = squidpy_p_full,
                          "Squidpy.means" = squidpy_p_full,
                          "NATMI.EPW" = natmi_epw_full,
                          "NATMI.SW" = natmi_sw_full,
                          "iTALK" = italk_full,
                          "Connectome.EPW" = conn_epw_full,
                          "Connectome.SW" = conn_sw_full,
                          "SCA" = sca_full)) %>%
    get_rank_frequencies()

# 5. PCA by Rank Frequencies
plot_freq_pca(rank_frequencies)




# Read NES
bc_nes <- read.csv("~/Repos/ligrec_decoupleR/input/sc_bc/breast_cancer_NES.csv",
                   header = TRUE, row.names = 1) %>%
    `colnames<-`(str_glue("c{seq(0, 11)}")) %>%
    `rownames<-`(str_glue("c{seq(0, 11)}"))


bc_nes_vec <- map2(.x=str_glue("c{seq(0, 11)}"), .y=seq(1, 12), .f=function(x1, y1){
    map2(.x=str_glue("c{seq(0, 11)}"), .y=seq(1, 12), .f=function(x2, y2){
        str_glue("{x1}_{x2}:{bc_nes[y1, y2]}")
        })
    }) %>%
    enframe() %>%
    unnest(value) %>%
    unnest(value) %>%
    separate(value, into = c("clust_pair", "NES"), sep = ":") %>%
    distinct() %>%
    select(clust_pair, NES) %>%
    mutate(NES = as.double(NES))
# mutate(NES = scale(as.double(NES))[,1])



# Ranked Frequencies x NES
rank_nes_freq <- rank_frequencies %>%
    left_join(., bc_nes_vec, by = "clust_pair") # %>%
# separate(clust_pair, into=c("c1","c2")) %>%
# filter(c1!=c2) %>%
# unite(c1, c2, col="clust_pair")

#lm
rank_nes_regression <- rank_nes_freq %>%
    group_by(name) %>%
    do(model = lm(freq ~ NES, data = .)) %>%
    mutate(adjr =  model %>% glance() %>% pluck("adj.r.squared")) %>%
    mutate(pval =  model %>% glance() %>% pluck("p.value")) %>%
    select(name, adjr, pval)  %>%
    separate(name, into = c("Method", "Resource"), convert = TRUE, sep = "_") %>%
    mutate_if(is.character, as.factor)


ggplot(rank_nes_regression, aes(x=adjr, y=-log10(pval), colour = Method, shape = Resource)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nlevels(rank_nes_regression$Method))) +
    scale_shape_manual(values=1:nlevels(rank_nes_regression$Resource)) +
    xlab("Adj. Rsq") +
    ggtitle("Linear Regression of Activities x NES")


# pearson corr
rank_nes_corr <- rank_nes_freq %>%
    # filter(str_detect(name, "OmniPath")) %>%
    group_by(name) %>%
    do(corr = cor.test(x = .$freq, y = .$NES, method = "pearson")) %>%
    mutate(coef = corr %>% glance() %>% pull(estimate),
           pval = corr %>% glance() %>% pull(p.value)) %>%
    select(name, coef, pval)  %>%
    separate(name, into = c("Method", "Resource"), convert = TRUE, sep = "_") %>%
    mutate_if(is.character, as.factor)

ggplot(rank_nes_corr, aes(x=coef, y=-log10(pval),
                          colour = Method, shape = Resource)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nlevels(rank_nes_corr$Method))) +
    scale_shape_manual(values=1:nlevels(rank_nes_regression$Resource)) +
    xlab("Pearson Correlation Coefficient")  # +
# ggtitle("Correlation of Cell-Pair Activities x NES")




# Binarized Activities
bc_freq <- sig_list %>%
    get_binary_frequencies()


# join NES and activities
bc_nes_freq <- bc_freq %>%
    left_join(., bc_nes_vec, by = "clust_pair")


# regression of Activity x NES
bc_nes_regression = bc_nes_freq %>%
    separate(clust_pair, into = c("c1", "c2")) %>%
    filter(c1!=c2) %>%
    unite("c1", "c2", col = "clust_pair") %>%
    group_by(name) %>%
    do(model = lm(freq ~ NES, data = .)) %>%
    mutate(adjr =  model %>% glance() %>% pluck("adj.r.squared")) %>%
    mutate(pval =  model %>% glance() %>% pluck("p.value")) %>%
    select(name, adjr, pval)  %>%
    separate(name, into = c("Method", "Resource"), convert = TRUE) %>%
    mutate_if(is.character, as.factor)


ggplot(bc_nes_regression, aes(x=adjr, y=-log10(pval), colour = Method, shape = Resource)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=brewer.pal(6, "Dark2")) +
    scale_shape_manual(values=1:nlevels(bc_nes_regression$Resource)) +
    xlab("Adj. Rsq")
