# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat323.rds")

# Get clusters
clusts <- breast_cancer@meta.data %>%
    select(seurat_clusters) %>%
    rownames_to_column("barcode") %>%
    as_tibble()

# Get coordinates
coords <- breast_cancer@images$slice1@coordinates  %>%
    rownames_to_column("barcode") %>%
    as_tibble()

# Join clusters and Coordinates
clust_coords <- coords %>%
    left_join(clusts, by = "barcode") %>%
    # calculate centroids
    group_by(seurat_clusters) %>%
    mutate(centroid_x = mean(row),
           centroid_y = mean(col))

# Plot Cell Coordinates
ggplot(clust_coords, aes(x=row, y=col, color = seurat_clusters)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(clust_coords$seurat_clusters)))) +
    scale_shape_manual(values=1:nlevels(clust_coords$seurat_clusters))

# Plot Cluster Centroid Coordinates
ggplot(clust_coords, aes(x=centroid_x, y=centroid_y, color = seurat_clusters)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(clust_coords$seurat_clusters)))) +
    scale_shape_manual(values=1:nlevels(clust_coords$seurat_clusters))


# Get Cluster Centroid info
clust_centroids <- clust_coords %>%
    select(seurat_clusters, centroid_x, centroid_y) %>%
    ungroup() %>%
    distinct()
clust_centroids


# Get Euclidean Distance between Centroids of different Clusters
centroid_combs <- clust_centroids %>%
    mutate(centroid_x2 = centroid_x,
           centroid_y2 = centroid_y) %>%
    mutate(clust1 = seurat_clusters,
           clust2 = seurat_clusters) %>%
    # Get all cluster combinatons
    tidyr::expand(clust1, clust2) %>%
    # join Centroid coordinates
    left_join(clust_centroids, by=c("clust1"="seurat_clusters")) %>%
    left_join(clust_centroids, by=c("clust2"="seurat_clusters")) %>%
    select(clust1,
           clust2,
           x1 = centroid_x.x,
           y1 = centroid_y.x,
           x2 = centroid_x.y,
           y2 =centroid_y.y)

# Centroid Euclidean Distance
centroid_eucl <- centroid_combs %>%
    # filter(clust1!=clust2) %>% # filter same cells
    mutate(eucl = sqrt((x1 - x2)**2 + (y1 - y2)**2)) %>%
    mutate(norm_eucl = scale(eucl)[,1]*-1) %>%
    unite(clust1, clust2, col = "clust_pair") %>%
    select(clust_pair, eucl, norm_eucl)
centroid_eucl

# Euc X Rank Frequencies
rank_euc_freq <- rank_frequencies %>%
    left_join(., centroid_eucl, by = "clust_pair")

rank_euc_corr <- rank_euc_freq %>%
    group_by(name) %>%
    na.omit() %>%
    do(corr = cor.test(x = .$freq, y = .$norm_eucl, method = "spearman")) %>%
    mutate(coef = corr %>% glance() %>% pull(estimate),
           pval = corr %>% glance() %>% pull(p.value)) %>%
    select(name, coef, pval)  %>%
    separate(name, into = c("Method", "Resource"), convert = TRUE, sep = "_") %>%
    mutate_if(is.character, as.factor)

ggplot(rank_euc_corr, aes(x=coef, y=-log10(pval), colour = Method, shape = Resource)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=brewer.pal(8, "Dark2")) +
    scale_shape_manual(values=1:nlevels(rank_euc_corr$Resource)) +
    xlab("Spearman Correlation Coefficient")  # +
# ggtitle("Correlation of Cell-Pair Activities x NES")
