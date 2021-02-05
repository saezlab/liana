library(msigdbr)
library(tidyverse)

# Hallmark GOs from MSig -------------------------------------------------------
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

ggplot(plot_data, aes(fill=gs_name, y=n, x=name)) +
    geom_bar(position="stack", stat="identity")

# These go terms are far too general and in this case
# in order to fill many NAs, I look at GO terms true for either the
# source or target (which means that GOs are likely to be false for either,
# since they generally match only one).
# Looking at ligs and receptors seperately results in ~50% NAs, while I already
# look at ~50 GO terms... a bit pointless.


# Cell Signalling Specific GOs ------------------------------------------------
# msigdb <- import_omnipath_annotations(resource = 'MSigDB', wide = TRUE)
BiocManager::install("GOfuncR")

library("GOfuncR")

# Get Cell Communication GOs and manually filter
ccc_gos <- GOfuncR::get_child_nodes("GO:0007154") %>%
    filter(distance <= 3) %>%
    select(child_go_id, child_name, distance) %>%
    filter(!grepl("cell communication", child_name)) %>%
    filter(!grepl("pollen", child_name))


# This would certainly be a useful thing to look at but we need more context
# i.e. if we are looking at a specific case study, we could use
# child GO terms of (e.g. kidney development, immune response, etc)
# to see which resource might be most approriate.

# Look at hierarchy and manually pick the most approriate terms:
# http://www.informatics.jax.org/vocab/gene_ontology/GO:0007154





# Pathway annotations from OmniPath
netpath_pathways <- import_omnipath_annotations(resource = 'NetPath', wide = TRUE)
# msigdb <- import_omnipath_annotations(resource = 'MSigDB', wide = TRUE) %>% filter(collection == 'hallmark')
signor_pathways <- import_omnipath_annotations(resource = 'SIGNOR', wide = TRUE)
signalink_pathways <- import_omnipath_annotations(resource = 'SignaLink_pathway', wide = TRUE)
