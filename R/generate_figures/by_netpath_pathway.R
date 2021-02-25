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
# Import intercell network and, since we cannot be
# sure of annotations attributed to complexes, remove them.
omnipath_intercell <- OmnipathR::import_intercell_network() %>% 
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))

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

pdf(file = "./figures/netpath_interaction.pdf",  
    width = 7,
    height = 5) 

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (NetPath)")

dev.off()

# RECEPTOR:
add_omnipath <- intercell_pathways_receptor %>% group_by(pathway) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_pathways_receptor %>% tidyr::separate_rows(sources, sep = ";") %>%
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

pdf(file = "./figures/netpath_receptor.pdf",  
    width = 7,
    height = 5) 

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Receptor pathway (NetPath)")

dev.off()

