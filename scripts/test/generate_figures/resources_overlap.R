require(UpSetR)
require(grid)
require(purrr)
source("./R/generate_figures/support_functions.R")

ligrec_interaction_list <- readRDS("./R/networks.RData")
upset_input <- make_upset_list(ligrec_interaction_list)

# Upset plot for omnipath and all resources within it:
pdf(file = "./figures/interactions_upset.pdf",  
    width = 9,
    height = 6) 
upset(fromList(upset_input), order.by = "freq", 
      nsets = length(upset_input), nintersects=30,
      scale.intersections = "log10") # Plots intersections on a log10 scale
dev.off()

# Upset plot for all non-empty resources within omnipath, excluding omnipath:
upset(fromList(within(upset_input, rm('OmniPath'))), order.by = "freq", nsets = length(upset_input)-1, nintersects=30)

# # Interactions which have consensus as being activating
# activating_interactions <- lapply(ligrec_interaction_list, function(x){
#   if (typeof(x)=='list') dplyr::filter(x, consensus_stimulation == 1)
#   else return(integer(0))
# }) %>% make_upset_list
# upset(fromList(activating_interactions), order.by = "freq", nsets = length(upset_input), nintersects=30)
# 
# 
# # Interactions which have consensus as being inactivating
# inactivating_interactions <- lapply(ligrec_interaction_list, function(x){
#   if (typeof(x)=='list') dplyr::filter(x, consensus_inhibition == 1)
#   else return(integer(0))
# }) %>% make_upset_list
# upset(fromList(inactivating_interactions), order.by = "freq", nsets = length(upset_input), nintersects=30)


