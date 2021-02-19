# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by Netpath pathway.
# Only interactions whose ligand and receptor belong to the same pathway,
# and only the top 15-16 pathways by cardinality are displayed.
# We also make an analogous figure with interactions between complexes removed.
require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./R/generate_figures/support_functions.R")

omni_list <- get_resource_list()
omnipath_intercell <- OmnipathR::import_intercell_network()
netpath_pathways <- import_omnipath_annotations(resource = 'NetPath', wide = TRUE) %>%
  dplyr::select(uniprot, pathway)

# top_pathways <- as.data.frame(summary(as.factor(netpath_pathways$pathway))) %>% 
#   setNames("count") %>%
#   dplyr::filter(count > 200)
# 
# top_netpath_pathways <- netpath_pathways %>% dplyr::filter(pathway %in% rownames(top_pathways))

intercell_pathways_receptor <- left_join(omnipath_intercell, netpath_pathways, by = c("target" = "uniprot"))
intercell_pathways_ligand <- left_join(omnipath_intercell, netpath_pathways, by = c("source" = "uniprot"))
# Restrict to interactions whose ligand and receptor belong to the same pathway
intercell_pathways <- intersect(intercell_pathways_receptor, intercell_pathways_ligand)

add_omnipath <- intercell_pathways %>% group_by(pathway) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_pathways %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, pathway) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, pathway, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(pathway))
top_categories <- data_w_omnipath %>% 
  group_by(pathway) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath %>%
  dplyr::filter(pathway %in% top_categories$pathway)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (NetPath)")

# The same but without complexes:
intercell_rm_complex <- omnipath_intercell %>% 
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))
intercell_pathways_receptor <- left_join(intercell_rm_complex, netpath_pathways, by = c("target" = "uniprot"))
intercell_pathways_ligand <- left_join(intercell_rm_complex, netpath_pathways, by = c("source" = "uniprot"))
# Restrict to interactions whose ligand and receptor belong to the same pathway
intercell_pathways_rm_complex <- intersect(intercell_pathways_receptor, intercell_pathways_ligand)

add_omnipath_rm_complex <- intercell_pathways_rm_complex %>% group_by(pathway) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath_rm_complex <- intercell_pathways_rm_complex %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, pathway) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath_rm_complex)%>%
  tidyr::complete(sources, pathway, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(pathway))
top_categories <- data_w_omnipath_rm_complex %>% 
  group_by(pathway) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath_rm_complex %>%
  dplyr::filter(pathway %in% top_categories$pathway)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (SignaLink)")

