# I) Enrichment by Cell Pair
tmp <- rank_frequencies %>%
    filter(str_detect(name, "_OmniPath")) %>%
    separate(clust_pair, into=c("source", "target"), sep = "_") %>%
    separate(name, into=c("tool", "resource"), sep = "_")  %>%
    select(tool, source, target, freq)

tmp1 <- tmp %>% select(tool, cp=source, freq) %>% mutate(cat="source")
tmp2 <- tmp %>% select(tool, cp=target, freq) %>% mutate(cat="target")

xd <- bind_rows(tmp1, tmp2)

# II) Enrichment by Resource
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


