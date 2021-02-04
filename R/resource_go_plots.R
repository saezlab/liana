library(msigdbr)
library(tidyverse)

# Hallmark GOs from MSig
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
go_format <- h_gene_sets %>%
    select("symbol" = human_gene_symbol,
           gs_name) %>%
    mutate(gs_name = str_replace(gs_name,"HALLMARK_", ""))

### Get all resource summaries by GO
## Get source GOs
source_go_summary <- omni_resources %>%
    map(function(op_resource){
        op_resource %>%
            # filter complexes
            filter(entity_type_intercell_source != "complex",
                   entity_type_intercell_target != "complex",) %>%
            select("source" = source_genesymbol,
                   "target" = target_genesymbol) %>%
            unite("interaction_symbols", c(source, target), remove = FALSE) %>%
            distinct_at(vars(interaction_symbols), .keep_all = TRUE) %>%
            distinct_at(vars(source), .keep_all = TRUE) %>%
            left_join(go_format, by = c("source" = "symbol"))  %>%
            group_by(gs_name) %>%
            summarise(n=n())
    }) %>%
    setNames(names(omni_resources)) %>%
    enframe() %>%
    unnest(value)


## Get target GO summary
target_go_summary <- omni_resources %>%
    map(function(op_resource){
        op_resource %>%
            # filter complexes
            filter(entity_type_intercell_source != "complex",
                   entity_type_intercell_target != "complex",) %>%
            select("source" = source_genesymbol,
                   "target" = target_genesymbol) %>%
            unite("interaction_symbols", c(target, target), remove = FALSE) %>%
            distinct_at(vars(interaction_symbols), .keep_all = TRUE) %>%
            distinct_at(vars(target), .keep_all = TRUE) %>%
            left_join(go_format, by = c("target" = "symbol"))  %>%
            group_by(gs_name) %>%
            summarise(n=n())
    }) %>%
    setNames(names(omni_resources)) %>%
    enframe() %>%
    unnest(value)

# Get Interaction GO summaries
interaction_go_summary <- omni_resources %>%
    map(function(op_resource){
        op_resource %>%
            # filter complexes
            filter(entity_type_intercell_source != "complex",
                   entity_type_intercell_target != "complex",) %>%
            select("source" = source_genesymbol,
                   "target" = target_genesymbol) %>%
            unite("interaction_symbols", c(source, target), remove = FALSE) %>%
            distinct_at(vars(interaction_symbols), .keep_all = TRUE) %>%
            distinct_at(vars(source), .keep_all = TRUE) %>%
            left_join(go_format, by = c("source" = "symbol")) %>%
            left_join(go_format, by = c("target" = "symbol")) %>%
            mutate(across(everything(), ~replace_na(.x, ""))) %>%
            unite("gs_name", c(gs_name.x, gs_name.y), sep = ", ")  %>%
            select(interaction_symbols, gs_name) %>%
            # correct format
            mutate(gs_name = gsub("^, $", "NAN", gs_name)) %>%
            mutate(gs_name = gsub("^, ", "", gs_name)) %>%
            mutate(gs_name = gsub(", $", "", gs_name)) %>%
            separate_rows(gs_name, sep=", ") %>%
            # summarize
            group_by(gs_name) %>%
            summarise(n=n())
    }) %>%
    setNames(names(omni_resources)) %>%
    enframe() %>%
    unnest(value)


top_interactions <- interaction_go_summary %>%
    filter(name=="OmniPath") %>%
    arrange(desc(n)) %>%
    top_n(10) %>%
    pull(gs_name)


plot_data <- interaction_go_summary %>%
    filter(gs_name %in% top_interactions)

# xd <- interaction_go_summary %>%
#     pivot_wider(id_cols = gs_name, names_from = name, values_from = n) %>%
#     mutate_at(vars(everything()), ~ replace(., is.na(.), 0)) %>%
#     arrange(desc(OmniPath)) %>%
#     column_to_rownames("gs_name") %>%
#     filter(sum(.$omnipath) < 100) %>%
#     rownames_to_column("gs_name")

ggplot(plot_data, aes(fill=gs_name, y=n, x=name)) +
    geom_bar(position="stack", stat="identity")
