# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by Signalink pathway.
# Only interactions whose ligand and receptor belong to the same pathway are displayed.
# We also make an analogous figure with interactions between complexes removed.

require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./scripts/test/generate_figures/support_functions.R")

omni_list <- get_resource_list()
# Import intercell network and, since we cannot be
# sure of annotations attributed to complexes, remove them.
omnipath_intercell <- OmnipathR::import_intercell_network() %>%
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))

signalink_pathways <- import_omnipath_annotations(resource = 'SignaLink_pathway', wide = TRUE) %>%
  dplyr::select(uniprot, pathway)

intercell_pathways_receptor <- left_join(omnipath_intercell, signalink_pathways, by = c("target" = "uniprot"))
intercell_pathways_ligand <- left_join(omnipath_intercell, signalink_pathways, by = c("source" = "uniprot"))
# Only interactions whose ligand and receptor belong to the same Signalink pathway
intercell_pathways <- intersect(intercell_pathways_ligand, intercell_pathways_receptor)

add_omnipath <- intercell_pathways %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  dplyr::mutate(sources = "OmniPath")
data_w_omnipath <- intercell_pathways %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>%
  group_by(sources, pathway) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, pathway, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(pathway))
#No need to filter for top categories because there are only 10

pdf(file = "./figures/signalink_interaction.pdf",
    width = 7,
    height = 5)

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=pathway, y=n, x=sources)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction pathway (SignaLink)")

dev.off()

# RECEPTOR:
add_omnipath <- intercell_pathways_receptor %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  dplyr::mutate(sources = "OmniPath")
data_w_omnipath <- intercell_pathways_receptor %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>%
  group_by(sources, pathway) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath) %>%
  tidyr::complete(sources, pathway, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(pathway))
#No need to filter for top categories because there are only 10

pdf(file = "./figures/signalink_receptor.pdf",
    width = 7,
    height = 5)

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=pathway, y=n, x=sources)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  labs(x = "Resource", y = "Number of interactions", fill = "Receptor pathway (SignaLink)")

dev.off()


# Each category (lig, rec, interaction) combined -------------------------------
# OmniPath split by intercell resources
# receptors reformatted
recs <- intercell_pathways_receptor %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>%
  group_by(sources, pathway) %>%
  filter(!is.na(pathway)) %>%
  summarise(n=n()) %>%
  tidyr::complete(sources, pathway, fill = list(n = 0))  %>%
  dplyr::filter(!is.na(pathway)) %>%
  ungroup() %>%
  mutate(type = "Receptor")

# ligs reformatted
ligs <- intercell_pathways_ligand %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>%
  group_by(sources, pathway) %>%
  summarise(n=n()) %>%
  tidyr::complete(sources, pathway, fill = list(n = 0))  %>%
  dplyr::filter(!is.na(pathway)) %>%
  ungroup() %>%
  mutate(type = "Ligand")

# interactions
interacts <- intercell_pathways %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>%
  group_by(sources, pathway) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  tidyr::complete(sources, pathway, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(pathway))  %>%
  mutate(type = "Interaction")


# Composite OmniPath
# omni receptors
omni_receptor <- intercell_pathways_receptor %>%
  select(source, target, pathway) %>%
  distinct() %>%
  group_by(pathway) %>%
  summarise(n=n()) %>%
  filter(!is.na(pathway)) %>%
  mutate(sources = "OmniPath") %>%
  mutate(type = "Receptor")

# omni_ligands
omni_ligand <- intercell_pathways_ligand %>%
  select(source, target, pathway) %>%
  distinct() %>%
  group_by(pathway) %>%
  summarise(n=n()) %>%
  filter(!is.na(pathway)) %>%
  mutate(sources = "OmniPath") %>%
  mutate(type = "Ligand")

# omni_interactions
omni_interacts <- intercell_pathways %>%
  select(target, source, pathway) %>%
  distinct() %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  dplyr::mutate(sources = "OmniPath") %>%
  mutate(type = "Interaction") %>%
  dplyr::filter(!is.na(pathway))

bind_all <- bind_rows(recs,
                      ligs,
                      interacts,
                      omni_ligand,
                      omni_receptor,
                      omni_interacts
                      ) %>%
  distinct() %>%
  rename(counts = n)


signalink_bar <- ggplot(bind_all, aes(fill=pathway, y=counts, x=type)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "Number of interactions", fill = "Pathways (SignaLink)") +
  facet_wrap( ~ sources, nrow = 1., scales ="free_x") +
  theme(panel.spacing = unit(0.2, "lines"), text = element_text(size=14))
signalink_bar

