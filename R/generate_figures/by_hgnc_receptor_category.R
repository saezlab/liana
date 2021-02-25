# This script generates a stacked bar plot to look at the intercell reactions
# contained in each resource by HGNC receptor category.
# Only interactions in the top 15-16 receptor categories by cardinality are displayed.
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

hgnc_receptor_categories <- import_omnipath_intercell(parent = 'receptor', 
                                                      resources = 'HGNC', 
                                                      scope = 'specific')

intercell_receptor_category <- left_join(omnipath_intercell, hgnc_receptor_categories, 
                                         by = c("target" = "uniprot"))
add_omnipath <- intercell_receptor_category %>% 
  group_by(category) %>% 
  summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
data_w_omnipath <- intercell_receptor_category %>% 
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
  labs(x = "Resource", y = "Number of interactions", fill = "Receptor category (HGNC)")
