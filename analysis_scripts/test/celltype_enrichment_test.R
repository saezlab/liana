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

