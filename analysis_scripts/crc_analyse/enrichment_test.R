# I) Enriched Cell Pairs
top_frac <- top_lists$top_500 %>%
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
         color = colorRampPalette(c("darkslategray2", "violetred2"))(100),
         fontsize = 18,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         cutree_cols = 7,
         treeheight_col = 100
         )


cellfraq_heat$tree_row$order
cellfraq_heat$tree_row$labels


# Get proportions of single cells and proper names for clusts
crc_form <- readRDS("input/crc_data/crc_korean_mod.rds")
crc_meta <- crc_form@meta.data
rm(crc_form)

cellcount_fraq <- crc_meta %>%
    select(Cell_clusters, Cell_subtype) %>%
    group_by(Cell_subtype) %>%
    summarise(cell_occur = n()) %>%
    mutate(cell_fraq = cell_occur/sum(cell_occur), .keep = "unused") %>%
    pivot_wider(names_from = Cell_subtype, values_from = cell_fraq) %>%
    mutate(mr = "Cell.Counts")


cellcount_heat <- pheatmap(cellcount_fraq %>%
                               column_to_rownames("mr") %>%
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
cellcount_heat




#### Enrichment of ligand-receptor within same cell type
# II) Enrichment by Method
full_resource <- compile_ligrec(omni_variants = TRUE, lr_pipeline = FALSE)

lr <- full_resource %>% ligrec_overlap()


lr_csea <- lr %>% ligand_receptor_classes('CancerSEA', state)
# lr_signal <- lr %>% ligand_receptor_classes('SignaLink_pathway', pathway, NULL)
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



# Get top 500 and bind csea states
top500_resource_tool <- get_swapped_list(top_lists$top_500)
top500_csea <- top500_resource_tool$OmniPath %>%
    enframe(name = "method", value = "results") %>%
    mutate(results = results %>% map(function(res) res %>%
                                         select(source, target, ligand, receptor))) %>%
    unnest(results) %>%
    pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
    pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
    left_join(csea_syms) %>%
    na.omit()

# Get enrichment per cell type by tool
top500_csea_enrich <-
    top500_csea %>%
    group_by(cell_type)
gk <- group_keys(top500_csea_enrich)

omni_csea_enrich <- top500_csea_enrich %>%
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
             "violetred2"))(100),
         annotation_col = annotations_df,
         annotation_colors = mycolors,
         fontsize = 18,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 100
)



# IV) Enriched by tool in the same system
top500_csea <- top500_resource_tool$OmniPath %>%
    enframe(name = "method", value = "results") %>%
    mutate(results = results %>% map(function(res) res %>%
                                         select(source, target, ligand, receptor))) %>%
    unnest(results) %>%
    pivot_longer(cols = c(ligand, receptor), names_to = "cat", values_to = "genesymbol") %>%
    pivot_longer(cols = c(source, target), names_to = "cell", values_to = "cell_type") %>%
    left_join(csea_syms) %>%
    na.omit()


# Get each resoure by tool
top500_system <- top_lists$top_500 %>%
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


gk <- group_keys(top500_system)

# Get enrichment on a "system" level, i.e. enrichment per tool by resource
top500_sys_enrich <- top500_system %>%
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


sys_heat_data <- top500_sys_enrich %>%
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
                          color = colorRampPalette(c("darkslategray2", "violetred2"))(100),
                          fontsize = 18,
                          drop_levels = TRUE,
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          border_color = NA,
                          treeheight_row = 0,
                          treeheight_col = 100
)






# V) Enrichment by Interactions (use CellChat DB)
conn_syms <- full_resource$OmniPath$connections %>%
    select(source, target, source_genesymbol, target_genesymbol)

csea_conns <- lr_csea$connections %>%
    select(source, target, state) %>%
    distinct() %>%
    left_join(conn_syms, by = c("source", "target")) %>%
    distinct() %>%
    select(source_genesymbol, target_genesymbol, state) %>%
    unite(source_genesymbol, target_genesymbol, col = "interaction")







# XXX I) Enrichment by Cell Pair
tmp <- rank_frequencies %>%
    filter(str_detect(name, "_OmniPath")) %>%
    separate(clust_pair, into=c("source", "target"), sep = "_") %>%
    separate(name, into=c("tool", "resource"), sep = "_")  %>%
    select(tool, source, target, freq)

tmp1 <- tmp %>% select(tool, cp=source, freq) %>% mutate(cat="source")
tmp2 <- tmp %>% select(tool, cp=target, freq) %>% mutate(cat="target")

xd <- bind_rows(tmp1, tmp2)



sign <- top_lists$top_500$Squidpy$CellPhoneDB %>%
    filter(source=="Proliferating" | target=="Proliferating") %>%
    pivot_longer(cols = c(ligand, receptor), values_to = "gene") %>%
    pull(gene) %>%
    unique()
hyp_obj <- hypeR(sign, genesets, background=2522)
hyp_obj



library(hypeR)
genesets <- msigdb_gsets("Homo sapiens", "H")

signature <- c("IDH3B","DLST","PCK2","CS","PDHB","PCK1","PDHA1","LOC642502",
               "PDHA2","LOC283398","FH","SDHD","OGDH","SDHB","IDH3A","SDHC",
               "IDH2","IDH1","OGDHL","PC","SDHA","SUCLG1","SUCLA2","SUCLG2")

hyp_obj <- hypeR(signature, genesets, background=2522)
hyp_obj









# XXX III) Denes
saveRDS(ligrec_form, "output/test/ligrec_form.RDS")
ligrec_form <- readRDS("output/test/ligrec_form.RDS")
ligrec_form
enrich_res <- enrich2(ligrec_form$connections,
                      var1 = pathway,
                      var2 = resource)
enrich_res



#' Take Pathway column and check if the occurrences of the Resource are
#' enriched with this pathway
var1="pathway"
var2="resource"

var1_str <- quo_text(sym(var1))
var2_str <- quo_text(sym(var2))

t_f <- c('TRUE', 'FALSE')

test_method = 'fisher'
p_adj_method = "fdr"

data <- ligrec_form$connections

data %<>% select(!!var1, !!var2)

data %>%
    {map(names(.), function(x){unique(.[[x]])})} %>%
    setNames(names(data)) %>%
    cross_df() %>%
    filter(!is.na(!!var1) & !is.na(!!var2)) %>%
    group_by(pathway, resource) %>%
    group_modify(
        function(.x, .y){
            f1 <- factor(data[[var1_str]] == .y[[var1_str]][1], levels = t_f)
            print(f1)
            f2 <- factor(data[[var2_str]] == .y[[var2_str]][1], levels = t_f)
            print(f2)
            result <- fisher.test(f1, f2)
            odds_ratio <- result$estimate
            if(test_method == 'barnard'){
                sink('NUL')
                result <- exec(barnard.test, !!!table(f1, f2), ...)
                sink()
            }else if(test_method == 'barnard2'){
                param <-
                    list(...) %>%
                    merge.list(list(fixed = NA, method = 'boschloo'))
                result <- exec(BarnardTest, f1, f2, !!!param)
            }
            tibble(pval = last(result$p.value), odds_ratio = odds_ratio)
        }
    ) %>%
    ungroup() %>%
    mutate(
        padj = p.adjust(pval, method = p_adj_method),
        enrichment = ifelse(
            odds_ratio < 1,
            -1 / odds_ratio,
            odds_ratio
        ) %>% unname
    )






# XXX II) Enrichment by Resource
# tmp <- sig_list_resource$OmniPath %>%
#     filter(str_detect(name, "_OmniPath")) %>%
#     separate(clust_pair, into=c("source", "target"), sep = "_") %>%
#     separate(name, into=c("tool", "resource"), sep = "_")  %>%
#     select(tool, source, target, freq)


# Get reso urces
full_resource <- compile_ligrec(omni_variants = TRUE, lr_pipeline = FALSE)

lr <- full_resource %>% ligrec_overlap()


lr_csea <- lr %>% ligand_receptor_classes('CancerSEA', state)
lr_signal <- lr %>% ligand_receptor_classes('SignaLink_pathway', pathway, NULL)



# bind symbols
lr_csea$connections
lr_signal$connections





lr_msig <- lr %>% ligand_receptor_classes(
    'MSigDB',
    geneset,
    filter_annot = collection == 'hallmark',
    label_annot = function(x){str_to_title(str_sub(x, 10))}
)


conn_syms <- full_resource$OmniPath$connections %>%
    select(source, target, source_genesymbol, target_genesymbol)


msig_conn_geneset <- lr_msig$connections %>%
    select(source, target, geneset) %>%
    distinct() %>%
    left_join(conn_syms, by = c("source", "target")) %>%
    distinct() %>%
    select(source_genesymbol, target_genesymbol, geneset) %>%
    unite(source_genesymbol, target_genesymbol, col = "interaction")



lig_syms <- full_resource$OmniPath$connections %>%
    select(source, target_genesymbol)



msig_lig_geneset <- lr_msig$ligands %>%
    select(source=uniprot, geneset) %>%
    distinct() %>%
    left_join(up_syms, by = c("source")) %>%
    distinct() %>%
    select(source_genesymbol, target_genesymbol, geneset) %>%
    unite(source_genesymbol, target_genesymbol, col = "interaction")


