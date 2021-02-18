library(OmnipathR)
library(UpSetR)
source("./R/generate_figures/support_functions.R")
ligrec_interaction_list <- readRDS("./R/networks.RData")

omnipath <- ligrec_interaction_list$Omnipath

# Finding actual categories:

omnipath_annotations <- import_omnipath_annotations(resources = names(ligrec_interaction_list))
omnipath_annotations_wider <- omnipath_annotations %>% 
  tidyr::pivot_wider(names_from = label, values_from = value) %>%
  group_by(uniprot) %>%
  summarise_all(funs(toString(unique(na.omit(.))))) %>%
  mutate_all(na_if,"")


role <- dplyr::select(omnipath_annotations_wider, c(uniprot, role))
location <- dplyr::select(omnipath_annotations_wider, c(uniprot, location))
mainclass <- dplyr::select(omnipath_annotations_wider, c(uniprot, mainclass))
receptor_class <- dplyr::select(omnipath_annotations_wider, c(uniprot, receptor_class))
# ... to be completed

add_omnipath_rc <- omnipath %>% 
  left_join(receptor_class, by = c("target"="uniprot")) %>% 
  group_by(receptor_class) %>% summarise(n = n()) %>%
  dplyr::mutate(sources = "Omnipath")

data_receptor_class <- omnipath %>% 
  left_join(receptor_class, by = c("target"="uniprot")) %>%
  tidyr::separate_rows(sources, sep = ";") %>%
  dplyr::filter(sources %in% names(ligrec_interaction_list)) %>% 
  group_by(sources, receptor_class) %>%
  summarise(n = n()) %>%
  ungroup %>%
  tidyr::complete(sources, receptor_class, fill = list(n = 0))

data_receptor_class_omnipath <- rbind(data_receptor_class, add_omnipath_rc)

# Stacked plot
ggplot(data_receptor_class, aes(fill=receptor_class, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Receptor Class")

# Stacked plot
ggplot(data_receptor_class_omnipath, aes(fill=receptor_class, y=n, x=sources)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Resource", y = "Number of interactions", fill = "Receptor Class")
########

cellphonedb <- ligrec_interaction_list$CellPhoneDB %>% dplyr::select(rowid, category_intercell_source, 
                                                                     category_intercell_target, source, target,
                                                                     genesymbol_intercell_source, 
                                                                     genesymbol_intercell_target)

# Location of receptor
cellphonedb_receptor_location <- cellphonedb %>% left_join(location, by = c("target" = "uniprot"))
#[1] "plasma membrane"           "plasma membrane, membrane" NA                          "plasma membrane, secreted"

# Location of ligand
cellphonedb_ligand_location <- cellphonedb %>% left_join(location, by = c("source" = "uniprot"))
#[1] "secreted"                            "plasma membrane, secreted, both"     NA                                   
#[4] "secreted, ecm"                       "plasma membrane, secreted"           "secreted, other"                    
#[7] "secreted, membrane"                  "secreted, both"                      "plasma membrane, secreted, membrane"
#[10] "plasma membrane"                     "plasma membrane, membrane"           "secreted, both, membrane"
  
# Receptor class
receptor_class <- omnipath_annotations_wider %>% dplyr::select(uniprot, receptor_class)
cellphonedb_receptor_class <- cellphonedb %>% left_join(receptor_class, by = c("target" = "uniprot"))
#[1] "chemokine_receptor_ccr"       "receptor"                     "cytokine_receptor"           
#[4] "atipical_chemokine_receptor"  "chemokine_receptor_cxcr"      "tnf_receptor"                
#[7] "estrogen_receptor"            "kir"                          "klr"                         
#[10] "growth_factor_receptor"       "hormone_receptor"             "inhibitory_lilrs"            
#[13] "artipical_chemokine_receptor" "cytokine_receptor_il6_family" "chemokine_receptor_xc"       
#[16] "chemokine_receptor_cx3cr1"    "chemokine_receptor"    

# Ligand class
secreted_class <- omnipath_annotations_wider %>% dplyr::select(uniprot, secreted_class)
cellphonedb_ligand_class <- cellphonedb %>% left_join(secreted_class, by = c("source" = "uniprot"))
#[1] "cytokine"              "immune-related"        "hormone"               "secreted"             
#[5] "cellsignal_wnt"        "growthfactor"          "cytokine;growthfactor" "growthfactor;hormone" 
#[9] "cytokine;hormone"      "cytokine_like"

# Pathway
pathway <-  omnipath_annotations_wider %>% dplyr::select(uniprot, pathway)
cellphonedb_ligand_pathway <- cellphonedb %>% left_join(pathway, by = c("source" = "uniprot"))
cellphonedb_receptor_pathway <- cellphonedb %>% left_join(pathway, by = c("target" = "uniprot"))
join_cols <- c("rowid", "category_intercell_source", "category_intercell_target", 
               "source", "target", "genesymbol_intercell_source", "genesymbol_intercell_target")
cellphonedb_pathway <- inner_join(cellphonedb_ligand_pathway, cellphonedb_receptor_pathway, by = join_cols, suffix = c("_ligand", "_receptor"))

#test <- cellphonedb_pathway %>% mutate(common_pathway = 
#                                         intersect(unlist(strsplit(pathway_ligand, ", ")), 
#                                                   unlist(strsplit(pathway_receptor, ", "))))

