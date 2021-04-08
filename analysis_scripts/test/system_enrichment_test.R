# IV) Enriched by tool in the same system
# Get each resoure by tool
top_system <- top_lists$top_250 %>%
    map(function(db){
        db %>%
            enframe(name = "resource", value = "results") %>%
            mutate(results = results %>% map(function(res) res %>%
                                                 select(source, target, ligand, receptor))) %>%
            unnest(results) %>%
            pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
            pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
            left_join(csea_syms) %>%
            na.omit()
    }) %>%
    enframe(name = "method", value = "results_resource") %>%
    unnest(results_resource) %>%
    group_by(resource)


gk <- group_keys(top_system)

# Get enrichment on a "system" level, i.e. enrichment per tool by resource
top_sys_enrich <- top_system %>%
    group_split() %>%
    setNames(gk$resource) %>%
    map(function(db) db %>%
            enrich2(var1 = state,
                    var2 = method)
    ) %>%
    enframe(name = "resource", value = "enrich") %>%
    unnest(enrich) %>%
    unite(method, resource, col="mr") %>%
    mutate_all(~ replace(., is.infinite(.), 0))


sys_heat_data <- top_sys_enrich %>%
    pivot_wider(id_cols = mr,
                names_from = state,
                values_from = enrichment,
                values_fill = 0)



# annotation groups (sequential vectors as in heatmap_binary_list)
method_groups <- sys_heat_data %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(method)
resource_groups <- sys_heat_data %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(resource)

# data frame with column annotations.
# with a column for resources and a column for methods
annotations_df <- data.frame(Resource = resource_groups,
                             Method = method_groups)  %>%
    mutate(rn = sys_heat_data$mr) %>%
    column_to_rownames("rn")

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
names(mycolors$Resource) <- unique(resource_groups)
names(mycolors$Method) <- unique(method_groups)



sys_enrich_heat <- pheatmap(sys_heat_data %>%
                                column_to_rownames("mr") %>%
                                t(),
                            annotation_col = annotations_df,
                            annotation_colors = mycolors,
                            display_numbers = FALSE,
                            silent = FALSE,
                            show_colnames = FALSE,
                            color = colorRampPalette(c("darkslategray2", "violetred2"))(20),
                            fontsize = 17,
                            drop_levels = TRUE,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            border_color = NA,
                            treeheight_row = 0,
                            treeheight_col = 100
)





# IV) Enriched by tool in the same system (signalink)
# Get each resoure by tool
top_signal_system <- top_lists$top_250 %>%
    map(function(db){
        db %>%
            enframe(name = "resource", value = "results") %>%
            mutate(results = results %>% map(function(res) res %>%
                                                 select(source, target, ligand, receptor))) %>%
            unnest(results) %>%
            pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
            pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
            left_join(signal_syms) %>%
            na.omit()
    }) %>%
    enframe(name = "method", value = "results_resource") %>%
    unnest(results_resource) %>%
    group_by(resource)


gk <- group_keys(top_signal_system)

# Get enrichment on a "system" level, i.e. enrichment per tool by resource
topsignal_sys_enrich <- top_signal_system %>%
    group_split() %>%
    setNames(gk$resource) %>%
    map(function(db) db %>%
            enrich2(var1 = pathway,
                    var2 = method)
    ) %>%
    enframe(name = "resource", value = "enrich") %>%
    unnest(enrich) %>%
    unite(method, resource, col="mr") %>%
    mutate_all(~ replace(., is.infinite(.), 0))


sys_heat_data <- topsignal_sys_enrich %>%
    pivot_wider(id_cols = mr,
                names_from = pathway,
                values_from = enrichment,
                values_fill = 0)




# annotation groups (sequential vectors as in heatmap_binary_list)
method_groups <- sys_heat_data %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(method)
resource_groups <- sys_heat_data %>%
    separate(mr, into = c("method", "resource"), sep = "_") %>%
    pull(resource)

# data frame with column annotations.
# with a column for resources and a column for methods
annotations_df <- data.frame(Resource = resource_groups,
                             Method = method_groups)  %>%
    mutate(rn = sys_heat_data$mr) %>%
    column_to_rownames("rn")

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
names(mycolors$Resource) <- unique(resource_groups)
names(mycolors$Method) <- unique(method_groups)

sys_enrich_heat <- pheatmap(sys_heat_data %>%
                                column_to_rownames("mr") %>%
                                t(),
                            annotation_col = annotations_df,
                            annotation_colors = mycolors,
                            display_numbers = FALSE,
                            silent = FALSE,
                            show_colnames = FALSE,
                            color = colorRampPalette(c("darkslategray2", "violetred2"))(20),
                            fontsize = 17,
                            drop_levels = TRUE,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            border_color = NA,
                            treeheight_row = 0,
                            treeheight_col = 100)



# VI) Enrichment with each method individually against whole geneset DB



# VII) Enrichment by Interactions (use CellChat DB)
conn_syms <- full_resource$OmniPath$connections %>%
    select(source, target, source_genesymbol, target_genesymbol)

csea_conns <- lr_csea$connections %>%
    select(source, target, state) %>%
    distinct() %>%
    left_join(conn_syms, by = c("source", "target")) %>%
    distinct() %>%
    select(source_genesymbol, target_genesymbol, state) %>%
    unite(source_genesymbol, target_genesymbol, col = "interaction")
