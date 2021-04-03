library(vegan)

# distance with Ranked results (a bit meaningless)
rank_dist <- rank_frequencies %>%
    group_by(name) %>%
    mutate(freq = range1_0(avg_rank)) %>% # convert freq to 0 to 1 (1 is highest res)
    separate(name, into = c("method", "resource"), sep = "_") %>%
    filter(resource == "OmniPath") %>%
    ungroup() %>%
    select(method, clust_pair, freq) %>%
    pivot_wider(id_cols = method, names_from = clust_pair , values_from = freq, values_fill = 0)  %>%
    column_to_rownames("method") %>%
    vegdist(method = "jaccard")
rank_dist

# range standardization (from 1 (highest rank) to 0 (lowest rank))
range1_0 <- function(x, ...){1 - (x - min(x, ...)) / (max(x, ...) - min(x, ...))}

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



tmp <- top_lists_swapped$top_500$OmniPath %>%
    prepForUpset() %>%
    select(-interaction) %>%
    t() %>%
    vegdist(method = "bray", diag=TRUE, upper = TRUE, binary = TRUE)




library(GGally)
ggpairs(tmp)


require(visNetwork, quietly = TRUE)

nodes <- data.frame(id = 1:3)
edges <- data.frame(from = c(1,2), to = c(1,3))
visNetwork(nodes, edges, width = "100%")


library(rpart)

# Basic classification tree
res <- rpart(Species~., data=iris)
visTree(res, main = "Iris classification Tree", width = "100%")



nodes <- data.frame(id = 1:15, label = paste("Label", 1:15),
                    group = sample(LETTERS[1:3], 15, replace = TRUE))

edges <- data.frame(from = trunc(runif(15)*(15-1))+1,
                    to = trunc(runif(15)*(15-1))+1)

edges$color <- "black"
edges$label <- NULL
edges$value <- c(1:15)
edges_bis <- edges
edges_bis$color <- "red"
# edges_bis$label <- "second"

visNetwork(nodes, rbind(edges, edges_bis, edges)) %>%
    visEdges(smooth = list(enabled = T, type = 'dynamic')) %>%
    visPhysics(solver = "barnesHut",
               barnesHut = list(springConstant = 0.002))


