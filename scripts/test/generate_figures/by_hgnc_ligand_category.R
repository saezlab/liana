# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by HGNC ligand category.
# Only interactions in the top 15-16 ligand categories by cardinality are displayed.
# We also make an analogous figure with interactions between complexes removed.
require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./R/generate_figures/support_functions.R")

omni_list <- get_resource_list()
omnipath_intercell <- OmnipathR::import_intercell_network()

hgnc_ligand_categories <- import_omnipath_intercell(parent = 'ligand', 
                                                    resources = 'HGNC', 
                                                    scope = 'specific')

# top_categories <- as.data.frame(summary(as.factor(hgnc_ligand_categories$category))) %>% 
#   setNames("count") %>%
#   dplyr::slice_max(order_by = count, n = 15)
# 
# ligands_in_top_categories <- hgnc_ligand_categories %>% 
#   dplyr::filter(category %in% rownames(top_categories))
# 
intercell_ligand_category <- left_join(omnipath_intercell, 
                                       hgnc_ligand_categories, 
                                       by = c("source" = "uniprot"))

# One dataframe representing all of Omnipath intercell
add_omnipath <- intercell_ligand_category %>% 
  group_by(category) %>% 
  summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
# Create dataframe for the individual resources AND Omnipath
# Remove NAs - corresponds to interactions which don't have an HGNC ligand category
data_w_omnipath <- intercell_ligand_category %>% 
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, category) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, category, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(category))

top_categories <- data_w_omnipath %>% 
  group_by(category) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath %>%
  dplyr::filter(category %in% top_categories$category)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=category, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Ligand category (HGNC)")

# The same but without complexes:
intercell_rm_complex <- omnipath_intercell %>% 
  dplyr::filter(!str_detect(source, "COMPLEX") & !str_detect(target, "COMPLEX"))

intercell_ligand_category_rm_complex <- left_join(intercell_rm_complex, 
                                                  hgnc_ligand_categories, 
                                                  by = c("source" = "uniprot"))

add_omnipath_rm_complex <- intercell_ligand_category_rm_complex %>% 
  group_by(category) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath_rm_complex <- intercell_ligand_category_rm_complex %>% 
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, category) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath_rm_complex)%>%
  tidyr::complete(sources, category, fill = list(n = 0)) %>%
  dplyr::filter(!is.na(category))
top_categories <- data_w_omnipath_rm_complex %>% 
  group_by(category) %>%
  summarise(count = sum(n)) %>%
  dplyr::slice_max(order_by = count, n = 15)
interactions_in_top_categories <- data_w_omnipath_rm_complex %>%
  dplyr::filter(category %in% top_categories$category)

# Stacked bar plot
ggplot(interactions_in_top_categories, aes(fill=category, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", 
       fill = "Ligand category (HGNC)")