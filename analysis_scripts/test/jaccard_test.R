library(vegan)

# range standardization (from 1 (highest rank) to 0 (lowest rank))
range1_0 <- function(x, ...){1 - (x - min(x, ...)) / (max(x, ...) - min(x, ...))}

# distance with Ranked results (a bit meaningless)
rank_dist <- rank_frequencies %>%
    group_by(name) %>%
    mutate(freq = range1_0(avg_rank)) %>% # convert freq to 0 to 1 (1 is highest res)
    separate(name, into = c("method", "resource"), sep = "_") %>%
    filter(resource == "OmniPath") %>%
    ungroup() %>%
    select(method, clust_pair, freq) %>%
    pivot_wider(id_cols = method,
                names_from = clust_pair,
                values_from = freq,
                values_fill = 0)  %>%
    column_to_rownames("method") %>%
    vegdist(method = "bray")
rank_dist



# data(BCI)
# H <- diversity(BCI)
# simp <- diversity(BCI, "simpson")
#
# vare.dist <- vegdist(varespec)
# vare.dist <- vegdist(decostand(varespec, "norm"), "euclidean")


# swap resource-tools
top_lists_swapped <-
    top_lists %>%
    map(function(tl) get_swapped_list(tl)) %>%
    setNames(names(top_lists))



tmp <- top_lists_swapped$top_2000$OmniPath %>%
    prepForUpset() %>%
    select(-interaction) %>%
    t() %>%
    vegdist(method = "jaccard", diag=FALSE, upper = TRUE, binary = TRUE) %>%
    as.matrix()
xy <- t(combn(colnames(tmp), 2))
df_pairs <- data.frame(xy, bray=tmp[xy]) %>% as_tibble() %>%
    dplyr::rename(method1="X1", method2="X2") %>%
    distinct()

nodes_x <- df_pairs %>%
    pull(method1) %>%
    unique() %>%
    as_tibble() %>%
    mutate(id=row_number()) %>%
    dplyr::rename(group = "value") %>%
    mutate(label = group)


edges_x <- df_pairs %>%
    arrange(method1, method2) %>%
    mutate(color = "blue",
           method1 = factor(method1),
           method2 = factor(method2)) %>%
    group_by(method1) %>%
    mutate(from = group_indices()) %>%
    ungroup() %>%
    group_by(method2) %>%
    mutate(to = group_indices() + 1) %>%
    mutate(value = 1 - bray) %>%
    mutate(title = "top2000")



###
tmp <- top_lists_swapped$top_50$OmniPath %>%
    prepForUpset() %>%
    select(-interaction) %>%
    t() %>%
    vegdist(method = "jaccard", diag=FALSE, upper = TRUE, binary = TRUE) %>%
    as.matrix()
xy <- t(combn(colnames(tmp), 2))
df_pairs <- data.frame(xy, bray=tmp[xy]) %>% as_tibble() %>%
    dplyr::rename(method1="X1", method2="X2") %>%
    distinct()



nodes_x <- df_pairs %>%
    pull(method1) %>%
    unique() %>%
    as_tibble() %>%
    mutate(id=row_number()) %>%
    dplyr::rename(group = "value") %>%
    mutate(label = group)


edges_x1 <- df_pairs %>%
    arrange(method1, method2) %>%
    mutate(color = "red",
           method1 = factor(method1),
           method2 = factor(method2)) %>%
    group_by(method1) %>%
    mutate(from = group_indices()) %>%
    ungroup() %>%
    group_by(method2) %>%
    mutate(to = group_indices() + 1) %>%
    mutate(value = 1 - bray) %>%
    mutate(title = "top500")

###
library(visNetwork)
visNetwork(nodes_x, bind_rows(edges_x, edges_x1)) %>%
    visEdges(smooth = list(enabled = TRUE, type = 'dynamic')) %>%
    visPhysics(solver = "forceAtlas2Based",
               forceAtlas2Based = list("gravitationalConstant" = -200,
                                       "springConstant" = 0.02))  %>%
    visLegend()




library(GGally)
ggpairs(tmp)
