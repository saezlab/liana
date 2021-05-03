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
