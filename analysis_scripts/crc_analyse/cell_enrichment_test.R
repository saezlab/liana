# I) Enriched Cell Pairs
top_frac <- top_lists$top_250 %>%
    map(function(db){
        db %>%
        enframe(name = "resource", value = "results") %>%
            mutate(results = results %>% map(function(res) res %>%
                                                 select(source, target))) %>%
            unnest(results) %>%
            pivot_longer(cols = c(source, target), names_to = "cat", values_to = "cell") %>%
            group_by(resource, cell) %>%
            summarise(cell_occur = n()) %>%
            mutate(cell_fraq = cell_occur/sum(cell_occur)) %>%
            ungroup() %>%
            select(resource, cell, cell_fraq) %>%
            pivot_wider(id_cols = resource,
                        names_from = cell,
                        values_from = cell_fraq,
                        values_fill = 0)
    }) %>%
    enframe(name = "method", value = "results_resource") %>%
    unnest(results_resource) %>%
    unite(method, resource, col = "mr")  %>%
    mutate_all(~ replace(., is.na(.), 0)) %>%
    select(-c("Stromal cells", "Epithelial cells"))


# annotation groups (sequential vectors as in heatmap_binary_list)
method_groups <- top_frac %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(method)
resource_groups <- top_frac %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(resource)

# data frame with column annotations.
# with a column for resources and a column for methods
annotations_df <- data.frame(Resource = resource_groups,
                             Method = method_groups)  %>%
    mutate(rn = top_frac$mr) %>%
    column_to_rownames("rn")

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
names(mycolors$Resource) <- unique(resource_groups)
names(mycolors$Method) <- unique(method_groups)


cellfraq_heat <- pheatmap(top_frac %>%
                              column_to_rownames("mr") %>%
                              t(),
         annotation_col = annotations_df,
         annotation_colors = mycolors,
         display_numbers = FALSE,
         silent = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("darkslategray2",
                                    "violetred2"))(20),
         fontsize = 18,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 100
         )


# Get proportions of single cells and proper names for clusts
crc_form <- readRDS("input/crc_data/crc_korean_mod.rds")
crc_meta <- crc_form@meta.data
rm(crc_form)
gc()

cellcount_fraq <- crc_meta %>%
    select(Cell_clusters, Cell_subtype) %>%
    group_by(Cell_subtype) %>%
    summarise(cell_occur = n()) %>%
    mutate(cell_fraq = cell_occur/sum(cell_occur), .keep = "unused") %>%
    pivot_wider(names_from = Cell_subtype, values_from = cell_fraq) %>%
    mutate(mr = "Cell.Counts")


cellfraq_heat$tree_row$labels
cellfraq_heat$tree_row$order

label_order <- data.frame(lab = cellfraq_heat$tree_row$labels,
           ord = cellfraq_heat$tree_row$order) %>%
    arrange(ord)


cellcount_heat <- pheatmap(cellcount_fraq %>%
                               column_to_rownames("mr") %>%
                               # select(label_order$lab) %>%
                               t(),
                           display_numbers = FALSE,
                           silent = FALSE,
                           show_colnames = FALSE,
                           color = colorRampPalette(c("darkslategray2", "violetred2"))(100),
                           fontsize = 18,
                           drop_levels = TRUE,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           border_color = NA,
                           treeheight_row = 0,
                           cutree_cols = 7,
                           treeheight_col = 100
                           )


