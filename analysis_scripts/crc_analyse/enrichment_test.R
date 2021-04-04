my_pallette = colorRampPalette(c(
    # "gray15",
    "darkslategray2",
    "violetred2"
    ))(100)


# I) Enriched Cell Pairs
# i.e. cells with the most edges that either stem from them or to them
top1000_resource_tool <- get_swapped_list(top_lists$top_1000)

cell_fraq <- top1000_resource_tool$OmniPath %>%
    enframe(name = "method", value = "results") %>%
    mutate(results = results %>% map(function(res) res %>%
            select(source, target))) %>%
    unnest(results) %>%
    pivot_longer(cols = c(source, target), names_to = "cat", values_to = "cell") %>%
    group_by(method, cell) %>%
    summarise(cell_occur = n()) %>%
    mutate(cell_fraq = cell_occur/sum(cell_occur)) %>%
    ungroup() %>%
    select(method, cell, cell_fraq)


cell_fraq %>%
    pivot_wider(id_cols = cell,
                names_from =  method,
                values_from = cell_fraq,
                values_fill = 0) %>%
    column_to_rownames("cell") %>%
    pheatmap()


#
top_frac <- top_lists$top_1000 %>%
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
            pivot_wider(id_cols = resource ,
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


pheatmap(top_frac %>% column_to_rownames("mr") %>% t(),
         annotation_col = annotations_df,
         annotation_colors = mycolors,
         display_numbers = FALSE,
         silent = FALSE,
         show_colnames = FALSE,
         color = my_pallette,
         fontsize = 18,
         drop_levels = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         treeheight_row = 0,
         treeheight_col = 100
         )


# II) Enrichment by Method




# III) Enrichment by Rank (e.g. PROGENy)








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




# XXX II) Enrichment by Resource
# tmp <- sig_list_resource$OmniPath %>%
#     filter(str_detect(name, "_OmniPath")) %>%
#     separate(clust_pair, into=c("source", "target"), sep = "_") %>%
#     separate(name, into=c("tool", "resource"), sep = "_")  %>%
#     select(tool, source, target, freq)

tmp1 <- tmp %>% select(tool, cp=source, freq) %>% mutate(cat="source")
tmp2 <- tmp %>% select(tool, cp=target, freq) %>% mutate(cat="target")

xd <- bind_rows(tmp1, tmp2)


full_resource <- compile_ligrec(omni_variants = TRUE, lr_pipeline = FALSE)

lr <- full_resource %>% ligrec_overlap()




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










import_omnipath_intercell(
    resource = 'HGNC',
    scope = 'specific',
    parent = "ligand",
    entity_type = 'protein'
) %>%
    select(genesymbol, category)

import_omnipath_intercell(
    resource = 'HGNC',
    scope = 'specific',
    parent = "receptor",
    entity_type = 'protein'
) %>%
    select(genesymbol, category)




lr_csea <- lr %>% ligand_receptor_classes('CancerSEA', state)
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



msig_lig_geneset


