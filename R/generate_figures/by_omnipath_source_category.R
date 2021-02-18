# This script generates a stacked bar plot to show the intercell interactions of
# each resource by source category. Then we create a summary dataframe to
# demonstrate that many of the interactions shown in the stacked bar plot are
# duplicates which contribute to multiple categories.

require(OmnipathR)
require(dplyr)
require(stringr)
require(ggplot2)
source("./R/generate_figures/support_functions.R")

omni_list <- get_resource_list()
omnipath_intercell <- OmnipathR::import_intercell_network()

# One dataframe representing all of Omnipath intercell
add_omnipath <- omnipath_intercell %>% 
  group_by(category_intercell_source) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")
# Create dataframe for the individual resources AND Omnipath
data_w_omnipath <- omnipath_intercell %>% 
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, category_intercell_source) %>%
  summarise(n = n()) %>%
  ungroup %>%
  rbind(add_omnipath)%>%
  tidyr::complete(sources, category_intercell_source, fill = list(n = 0))

# Stacked bar plot
ggplot(data_w_omnipath, aes(fill=category_intercell_source, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Source category")

# Exploring the number of duplicates belonging to each category:
data <- omnipath_intercell %>% tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, category_intercell_source) %>%
  summarise(n = n()) %>%
  ungroup %>%
  tidyr::complete(sources, category_intercell_source, fill = list(n = 0))

sum_data <- data %>% group_by(sources) %>%
  summarise(count = sum(n))

count_duplicates <- omnipath_intercell %>% 
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% omni_list) %>% 
  group_by(sources, source, target) %>%
  summarise(duplicates = n()) %>%
  ungroup %>%
  group_by(sources) %>%
  summarise(duplicated_interactions = sum(duplicates) - n(), 
            unique_interactions = n())

omnipath_duplicates <- omnipath_intercell %>% 
  group_by(source, target) %>% 
  summarise(duplicates = n()) %>%
  ungroup %>%
  summarise(duplicated_interactions = sum(duplicates) - n(), 
            unique_interactions = n()) %>%
  mutate(sources = "Omnipath", count = nrow(omnipath_intercell))

summary <- full_join(sum_data, count_duplicates) %>% 
  rbind(omnipath_duplicates) %>% 
  mutate(fraction_unique = unique_interactions/count)




