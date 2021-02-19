# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by Signalink pathway.
# Only interactions whose ligand and receptor belong to the same pathway are displayed.
# We also make an analogous figure with interactions between complexes removed.

require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./R/generate_figures/support_functions.R")

omni_list <- get_resource_list()
omnipath_intercell <- OmnipathR::import_intercell_network()

signalink_pathways <- import_omnipath_annotations(resource = 'SignaLink_pathway', wide = TRUE) %>%
  dplyr::select(uniprot, pathway)
  
intercell_pathways_receptor <- left_join(omnipath_intercell, signalink_pathways, by = c("target" = "uniprot"))
intercell_pathways_ligand <- left_join(omnipath_intercell, signalink_pathways, by = c("source" = "uniprot"))
# Only interactions whose ligand and receptor belong to the same Signalink pathway
intercell_pathways <- intersect(intercell_pathways_ligand, intercell_pathways_receptor)

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
#No need to filter for top categories because there are only 10

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (SignaLink)")

# The same but without complexes:
intercell_rm_complex <- omnipath_intercell %>% 
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))

intercell_pathways_receptor <- left_join(intercell_rm_complex, signalink_pathways, by = c("target" = "uniprot"))
intercell_pathways_ligand <- left_join(intercell_rm_complex, signalink_pathways, by = c("source" = "uniprot"))
# Only interactions whose ligand and receptor belong to the same Signalink pathway
intercell_pathways_rm_complex <- intersect(intercell_pathways_ligand, intercell_pathways_receptor)

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
#No need to filter for top categories because there are only 10
# Stacked bar plot
ggplot(data_w_omnipath_rm_complex, aes(fill=pathway, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (SignaLink)")

