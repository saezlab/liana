library(UpSetR)
ligrec_interaction_list <- readRDS("./networks.RData")

make_upset_list <- function(ligrec_list){
  lapply(ligrec_list, function(x){
  if(typeof(x) == "list"){
    unlist(dplyr::select(x,rowid), use.names=FALSE)
  }
  else{
    return(integer(0))
  }
  })
}
upset_input <- make_upset_list(ligrec_interaction_list)

# Upset plot for omnipath and all resources within it:
upset(fromList(upset_input), order.by = "freq", nsets = length(upset_input), nintersects=30)

# Upset plot for all non-empty resources within omnipath, excluding omnipath:
upset(fromList(within(upset_input, rm('Omnipath'))), order.by = "freq", nsets = length(upset_input)-1, nintersects=30)

# Interactions which have consensus as being activating
activating_interactions <- lapply(ligrec_interaction_list, function(x){
  if (typeof(x)=='list') dplyr::filter(x, consensus_stimulation == 1)
  else return(integer(0))
}) %>% make_upset_list
upset(fromList(activating_interactions), order.by = "freq", nsets = length(upset_input), nintersects=30)


# Interactions which have consensus as being inactivating
inactivating_interactions <- lapply(ligrec_interaction_list, function(x){
  if (typeof(x)=='list') dplyr::filter(x, consensus_inhibition == 1)
  else return(integer(0))
}) %>% make_upset_list
upset(fromList(inactivating_interactions), order.by = "freq", nsets = length(upset_input), nintersects=30)


