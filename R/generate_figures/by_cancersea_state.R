# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by CancerSEA state.
# Only interactions whose ligand and receptor belong to the same state are displayed.
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

cancersea_states <- import_omnipath_annotations(resource = 'CancerSEA', wide = TRUE)%>%
  dplyr::select(uniprot, state)

intercell_states_receptor <- left_join(omnipath_intercell, cancersea_states, by = c("target" = "uniprot"))
intercell_states_ligand <- left_join(omnipath_intercell, cancersea_states, by = c("source" = "uniprot"))
# Only interactions whose ligand and receptor belong to the same Signalink state
intercell_states <- intersect(intercell_states_ligand, intercell_states_receptor)

add_omnipath <- intercell_states %>% group_by(state) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_states %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, state) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, state, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(state))
#No need to filter for top categories because there are only 10

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=state, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Interaction state (CancerSEA)")

# Receptors (working on writing a function so code doesn't need to be duplicated)
add_omnipath <- intercell_states_receptor %>% group_by(state) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_states_receptor %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, state) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, state, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(state))
#No need to filter for top categories because there are only 10

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=state, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interaction", fill = "Receptor state (CancerSEA)")

