library(UpSetR)
library(grid)
source("./R/generate_figures/support_functions.R")

ligrec_interaction_list <- readRDS("./R/networks.RData")

# Doing this by uniprot ID, but it could also be done by gene ID
ligands_omnipath_uniprot <- unique(ligrec_interaction_list$Omnipath$source) %>% 
  as.data.frame %>% 
  setNames("uniprot") %>%
  tibble::rowid_to_column() 

ligands_list <- names(ligrec_interaction_list) %>% 
  map(function(x){
    message(x)
    if(x!="Omnipath" & typeof(ligrec_interaction_list[[x]]) == "list"){
      ligands <- unique(ligrec_interaction_list[[x]]$source) %>% as.data.frame %>% setNames("uniprot") %>%
        left_join(ligands_omnipath_uniprot, by = "uniprot")
    } else if(x == "Omnipath"){
      ligands_omnipath_uniprot
    } else {
      integer(0)
    }
  }) %>% setNames(names(ligrec_interaction_list))

ligands_input <- make_upset_list(ligands_list)
upset(fromList(ligands_input), order.by = "freq", nsets = length(ligands_input), nintersects=30)
grid.text("Ligands by Uniprot ID",x = 0.65, y=0.95, gp=gpar(fontsize=10))

# Same for receptors

receptors_omnipath_uniprot <- unique(ligrec_interaction_list$Omnipath$target) %>% 
  as.data.frame %>%
  setNames("uniprot") %>%
  tibble::rowid_to_column()

receptors_list <- names(ligrec_interaction_list) %>% 
  map(function(x){
    message(x)
    if(x!="Omnipath" & typeof(ligrec_interaction_list[[x]]) == "list"){
      receptors <- unique(ligrec_interaction_list[[x]]$target) %>% as.data.frame %>% setNames("uniprot") %>%
        left_join(receptors_omnipath_uniprot, by = "uniprot")
    } else if(x == "Omnipath"){
      receptors_omnipath_uniprot
    } else {
      integer(0)
    }
  }) %>% setNames(names(ligrec_interaction_list))

receptors_input <- make_upset_list(receptors_list)
upset(fromList(receptors_input), order.by = "freq", nsets = length(receptors_input), nintersects=30)
grid.text("Receptors by Uniprot ID",x = 0.65, y=0.95, gp=gpar(fontsize=10))

