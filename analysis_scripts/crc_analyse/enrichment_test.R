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


#### Enrichment of ligand-receptor within same cell type
# II) Enrichment by Method
full_resource <- compile_ligrec(omni_variants = TRUE, lr_pipeline = FALSE)

lr <- full_resource %>% ligrec_overlap()


lr_csea <- lr %>% ligand_receptor_classes('CancerSEA', state)
lr_signal <- lr %>% ligand_receptor_classes('SignaLink_pathway', pathway, NULL)
# lr_msig <- lr %>% ligand_receptor_classes(
#     'MSigDB',
#     geneset,
#     filter_annot = collection == 'hallmark',
#     label_annot = function(x){str_to_title(str_sub(x, 10))}
# )


# bind symbols
ligand_syms <- full_resource$OmniPath$ligands %>%
    select(uniprot, genesymbol)

receptor_syms <- full_resource$OmniPath$receptors %>%
    select(uniprot, genesymbol)


# CancerSEA
csea_lig_syms <- lr_csea$ligands %>%
    select(uniprot, state) %>%
    na.omit() %>%
    left_join(ligand_syms, by = "uniprot") %>%
    select(genesymbol, state) %>%
    distinct()

csea_rec_syms <- lr_csea$receptors %>%
    select(uniprot, state) %>%
    na.omit() %>%
    left_join(receptor_syms, by = "uniprot")  %>%
    select(genesymbol, state) %>%
    distinct()

csea_syms <- bind_rows(csea_lig_syms, csea_rec_syms)



# Signalink
signal_lig_syms <- lr_signal$ligands %>%
    select(uniprot, pathway) %>%
    na.omit() %>%
    left_join(ligand_syms, by = "uniprot") %>%
    select(genesymbol, pathway) %>%
    distinct()

signal_rec_syms <- lr_signal$receptors %>%
    select(uniprot, pathway) %>%
    na.omit() %>%
    left_join(receptor_syms, by = "uniprot")  %>%
    select(genesymbol, pathway) %>%
    distinct()

signal_syms <- bind_rows(signal_lig_syms, signal_rec_syms)


# Get swapped list
top250_resource_tool <- get_swapped_list(top_lists$top_250)


# Get top 250 and bind csea states
top250_csea <- top250_resource_tool$OmniPath %>%
    enframe(name = "method", value = "results") %>%
    mutate(results = results %>% map(function(res) res %>%
                                         select(source, target, ligand, receptor))) %>%
    unnest(results) %>%
    pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
    pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
    left_join(csea_syms) %>%
    na.omit()

# Get enrichment per cell type by tool
top250_csea_enrich <-
    top250_csea %>%
    group_by(cell_type)
gk <- group_keys(top250_csea_enrich)

omni_csea_enrich <- top250_csea_enrich %>%
    group_split() %>%
    map(function(met)
        met %>%
            enrich2(var1 = state,
                    var2 = method)) %>%
    setNames(gk$cell_type) %>%
    enframe(name = "cell_type", value = "results") %>%
    unnest(results) %>%
    filter(abs(enrichment) != Inf)

oce_heat <- omni_csea_enrich %>%
    select(method, cell_type, state, enrichment) %>%
    unite(method, cell_type, col = "mc", sep = "_") %>%
    pivot_wider(id_cols = mc, names_from = state, values_from = enrichment) %>%
    mutate_all(~ replace(., is.na(.), 0))


# annotation groups (sequential vectors as in heatmap_binary_list)
method_groups <- oce_heat %>%
    separate(mc, into = c("method", "cell"), sep = "_") %>%
    pull(method)
cell_types <- oce_heat %>%
    separate(mc, into = c("method", "cell"), sep = "_") %>%
    pull(cell)

# data frame with column annotations.
# with a column for resources and a column for methods
annotations_df <- data.frame(Cell_type = cell_types,
                             Method = method_groups)  %>%
    mutate(rn = oce_heat$mc) %>%
    column_to_rownames("rn")

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Cell_type = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(cell_types))))
names(mycolors$Cell_type) <- unique(cell_types)
names(mycolors$Method) <- unique(method_groups)

# heatmap
pheatmap(oce_heat %>%
             column_to_rownames("mc") %>%
             t(),
         display_numbers = FALSE,
         silent = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c(
             "darkslategray2",
             "violetred2"))(20),
         annotation_col = annotations_df,
         annotation_colors = mycolors,
         fontsize = 16,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 100
)




# SIGNALINK
top250_signal <- top250_resource_tool$OmniPath %>%
    enframe(name = "method", value = "results") %>%
    mutate(results = results %>% map(function(res) res %>%
                                         select(source, target, ligand, receptor))) %>%
    unnest(results) %>%
    pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
    pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
    left_join(signal_syms) %>%
    na.omit()

# Get enrichment per cell type by tool
top250_signal_enrich <-
    top250_signal %>%
    group_by(cell_type)
gk <- group_keys(top250_signal_enrich)

omni_csea_enrich <- top250_signal_enrich %>%
    group_split() %>%
    map(function(met)
        met %>%
            enrich2(var1 = pathway,
                    var2 = method)) %>%
    setNames(gk$cell_type) %>%
    enframe(name = "cell_type", value = "results") %>%
    unnest(results) %>%
    filter(abs(enrichment) != Inf)

oce_heat <- omni_csea_enrich %>%
    select(method, cell_type, pathway, enrichment) %>%
    unite(method, cell_type, col = "mc", sep = "_") %>%
    pivot_wider(id_cols = mc, names_from = pathway, values_from = enrichment) %>%
    mutate_all(~ replace(., is.na(.), 0))


# annotation groups (sequential vectors as in heatmap_binary_list)
method_groups <- oce_heat %>%
    separate(mc, into = c("method", "cell"), sep = "_") %>%
    pull(method)
cell_types <- oce_heat %>%
    separate(mc, into = c("method", "cell"), sep = "_") %>%
    pull(cell)

# data frame with column annotations.
# with a column for resources and a column for methods
annotations_df <- data.frame(Cell_type = cell_types,
                             Method = method_groups)  %>%
    mutate(rn = oce_heat$mc) %>%
    column_to_rownames("rn")

# List with colors for each annotation.
mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                 Cell_type = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(cell_types))))
names(mycolors$Cell_type) <- unique(cell_types)
names(mycolors$Method) <- unique(method_groups)

# heatmap
pheatmap(oce_heat %>%
             column_to_rownames("mc") %>%
             t(),
         display_numbers = FALSE,
         silent = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c(
             "darkslategray2",
             "violetred2"))(20),
         annotation_col = annotations_df,
         annotation_colors = mycolors,
         fontsize = 17,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 100
)






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



