# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by MSigDB Hallmark geneset.
# Only interactions whose ligand and receptor belong to the same geneset,
# and only the top 15-16 genesets by cardinality are displayed.
# We also make an analogous figure with interactions between complexes removed.
require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./R/generate_figures/support_functions.R")

omni_list <- get_resource_list()
omnipath_intercell <- OmnipathR::import_intercell_network()

msigdb <- import_omnipath_annotations(resource = 'MSigDB', wide = TRUE)

msigdb_hallmark <- msigdb %>% dplyr::filter(collection == 'hallmark') %>% 
  dplyr::mutate(geneset = str_replace(geneset, "HALLMARK_", "")) %>%
  dplyr::select(uniprot, geneset)

intercell_hallmark_geneset_receptor <- left_join(omnipath_intercell, msigdb_hallmark, by = c("target" = "uniprot"))
intercell_hallmark_geneset_ligand <- left_join(omnipath_intercell, msigdb_hallmark, by = c("source" = "uniprot"))
# Interactions whose ligand and receptor belong to the same Hallmark geneset:
intercell_hallmark_geneset_interactions <- intersect(intercell_hallmark_geneset_ligand, intercell_hallmark_geneset_receptor)

add_omnipath <- intercell_hallmark_geneset_interactions %>% group_by(geneset) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_hallmark_geneset_interactions %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, geneset) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, geneset, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(geneset))
top_categories <- data_w_omnipath %>% 
  group_by(geneset) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath %>%
  dplyr::filter(geneset %in% top_categories$geneset)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=geneset, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction geneset (MSigDB, Hallmark)")

# The same but without complexes:
intercell_rm_complex <- omnipath_intercell %>% 
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))

intercell_hallmark_geneset_receptor <- left_join(intercell_rm_complex, 
                                                 msigdb_hallmark, 
                                                 by = c("target" = "uniprot"))
intercell_hallmark_geneset_ligand <- left_join(intercell_rm_complex, 
                                               msigdb_hallmark, 
                                               by = c("source" = "uniprot"))
# Interactions whose ligand and receptor belong to the same Hallmark geneset:
intercell_hallmark_geneset_rm_complex <- intersect(intercell_hallmark_geneset_ligand, 
                                                     intercell_hallmark_geneset_receptor)

add_omnipath_rm_complex <- intercell_hallmark_geneset_rm_complex %>% 
  group_by(geneset) %>% 
  summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath_rm_complex <- intercell_hallmark_geneset_rm_complex %>% 
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, geneset) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath_rm_complex)%>%
  tidyr::complete(sources, geneset, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(geneset))
top_categories <- data_w_omnipath_rm_complex %>% 
  group_by(geneset) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath_rm_complex %>%
  dplyr::filter(geneset %in% top_categories$geneset)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=geneset, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction geneset (MSigDB, Hallmark)")
