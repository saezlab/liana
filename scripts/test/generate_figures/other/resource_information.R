library(OmnipathR)

annotations <- import_omnipath_annotations(resources = names(ligrec_interaction_list), wide = TRUE)
information_in_each_resource <- lapply(names(annotations), function(x){
  message(x)
  message(colnames(annotations[[x]]))
  return(colnames(annotations[[x]]))
}) %>% setNames(names(annotations))

library(qdapTools)
annotations_in_resource <- mtabulate(information_in_each_resource) %>% as.data.frame()
df <- annotations_in_resource %>%
  tibble::rownames_to_column() %>%
  tidyr::gather(colname, value, -rowname)

ggplot(df, aes(x = colname, y = rowname, fill = value)) +
  geom_tile()  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  labs(x = "Annotation", y = "Resource")

# Possibilities for categories:
adhesome_mainclasses <- levels(as.factor(annotations$Adhesome$mainclass))
#[1] "Actin regulation"             "Adaptor"                      "Adhesion receptor"           
#[4] "cAMP phosphodiesterase"       "Channel"                      "Chaperone"                   
#[7] "E3-ligase"                    "GAP"                          "GEF"                         
#[10] "GTPase"                       "Phospholipase"                "Protease"                    
#[13] "PtdIns kinase"                "PtdIns phosphatase"           "RNA or DNA regulation"       
#[16] "Serine palmitoyltransferase"  "Serine phosphatase"           "Serine/threonine kinase"     
#[19] "Serine/threonine phosphatase" "Tyrosine kinase"              "Tyrosine Kinase"             
#[22] "Tyrosine phosphatase"         "Unknown" 
baccin2019_subclasses <- levels(as.factor(annotations$Baccin2019$subclass))
#[1] "chemokine"              "chemokine_receptor"     "cytokine"               "cytokine_receptor"     
#[5] "growth_factor"          "growth_factor_receptor" "interleukin"            "interleukin_receptor"  
#[9] "other"  
cellchatdb_categories <- levels(as.factor(annotations$CellChatDB$category))
#[1] "Cell-Cell Contact"  "ECM-Receptor"       "Secreted Signaling"
cellphonedb_receptor_classes <- levels(as.factor(annotations$CellPhoneDB$receptor_class))
#[1] ""                             "artipical_chemokine_receptor" "atipical_chemokine_receptor" 
#[4] "chemokine_receptor"           "chemokine_receptor_ccr"       "chemokine_receptor_cx3cr1"   
#[7] "chemokine_receptor_cxcr"      "chemokine_receptor_xc"        "cytokine_receptor"           
#[10] "cytokine_receptor_il6_family" "estrogen_receptor"            "growth_factor_receptor"      
#[13] "hla"                          "hormone_receptor"             "inflammation"                
#[16] "inhibitory_lilrs"             "interferon_receptor"          "kir"                         
#[19] "klr"                          "receptor"                     "tgfbeta_receptor"            
#[22] "tnf_receptor"  
cellphonedb_secreted_classes <- levels(as.factor(annotations$CellPhoneDB$secreted_class))
#[1] ""                      "cellsignal_wnt"        "cytokine"              "cytokine_like"        
#[5] "cytokine;growthfactor" "cytokine;hormone"      "growthfactor"          "growthfactor;hormone" 
#[9] "hormone"               "immune-related"        "secreted"    
hmpr_mainclasses <- levels(as.factor(annotations$HPMR$mainclass))
#[1] "Cytokine Type 1 receptors"   "Cytokine Type 2 receptors"    "GPI-anchored"                                                            
#[4] "Guanylyl Cyclase receptors"  "Integrins"                    "Interleukin-17 receptors"                                                
#[7] "LINGO coreceptors for Nogo/p75"                                          
#[8] "Low-density lipoprotein (LDL) receptor and LDL receptor-related proteins"
#[9] "LRR-Ig Receptors"            "Netrin receptors"            "Neurexins"   
#[12] "Notch"                      "Other receptor/family"       "Patched"   
#[15] "Plexins"                    "Receptor Tyrosine Kinases (RTK)"                                         
#[17] "Receptor-like protein tyrosine phosphatases (RPTPs)"                     
#[18] "Roundabout"                 "Seven transmembrane (7TM) receptors"    
#[20] "Tetraspanins"               "TNF/NGF"                     "Toll"        
hmpr_subclasses <- levels(as.factor(annotations$HPMR$subclass))
#[1] "7TM A"                                      "7TM B"                                     
#[3] "7TM C"                                      "7TM D"                                     
#[5] "CD63"                                       "CICYTR (Cytokine Type 1 receptors)"        
#[7] "DCC receptors"                              "Discoidin domain receptor"                 
#[9] "EPHRIN"                                     "ERBB/EGF"                                  
#11] "FGR"                                        "Folate Receptors"                          
#[13] "Frizzled & Smoothened"                      "GHR (Cytokine Type 1 receptors)"           
#[15] "IL-1R2 (Cytokine Type 1 receptors)"         "IL-21R (Cytokine Type 1 receptors)"        
#[17] "IL2RA"                                      "IL2RG(Cytokine Type 1 receptors)"          
#[19] "Immune Cell Receptors"                      "INSULIN-R"                                 
#[21] "ITAB"                                       "ITAM"                                      
#[23] "ITB-4"                                      "Leukocyte Ig-like Receptors"               
#[25] "MET"                                        "Miscellaneous"                             
#27] "Neuropilins"                                "NGFR/NTR/TRK"                              
#[29] "Non-integrin laminin-binding proteins"      "OSMR(Cytokine Type 1 receptors)"           
#[31] "Phagocytosis receptors"                     "Prolactin Receptors"                       
#[33] "RAMP"                                       "ROR"                                       
#[35] "RPTPETA"                                    "RPTPOIC"                                   
#[37] "RPTPT12"                                    "RPTPUI5"                                   
#[39] "RPTPZETA"                                   "Scavenger"                                 
#[41] "Selectin Receptors"                         "Seven TM - Other"                          
#[43] "Subfamily1 (Netrin receptors)"              "Subfamily2 (Netrin receptors)"             
#[45] "Syndecan proteoglycan"                      "T13C (TNF/NGF)"                            
#[47] "TET3"                                       "TGF-beta serine/threonine kinase receptors"
#[49] "TIL (Toll)"                                 "TLR9 (Toll)"                               
#[51] "TM4SF5"                                     "TNR4 (TNF/NGF)"                            
#[53] "TNR8 (TNF/NGF)"                             "TNR9 (TNF/NGF)"                            
#[55] "TR12 (TNF/NGF)"                             "TR16 (TNF/NGF)"                            
#[57] "TRAIL (TNF/NGF)"                            "Transferrin Receptors"                     
#[59] "TYR/MER/UFO"                                "VEGF/PDGF" 
hmpr_subsubclasses <- levels(as.factor(annotations$HPMR$subsubclass))
icellnet_families <- levels(as.factor(annotations$ICELLNET$family))
#[1] "Antigen binding"  "Checkpoint"       "Chemokine"        "Cytokine"         "Growth factor"   
#[6] "Hormone"          "Neuropeptide"     "Notch signalling"
icellnet_classifications <- levels(as.factor(annotations$ICELLNET$classification))
#[1] "Antigen binding"                                 "Checkpoint"                                     
#[3] "Checkpoint;Cytokine"                             "Checkpoint;Cytokine;Cytokine type 1;Interleukin"
#[5] "Checkpoint;Cytokine;Tgf"                         "Checkpoint;Cytokine;Tnf family"                 
#[7] "Chemokine"                                       "Chemokine;Interleukin"                          
#[9] "Cytokine"                                        "Cytokine;Cytokine type 1"                       
#[11] "Cytokine;Cytokine type 1;Interleukin"            "Cytokine;Cytokine type 2"                       
#[13] "Cytokine;Cytokine type 2;Interleukin"            "Cytokine;Il1;Interleukin"                       
#[15] "Cytokine;Il17;Interleukin"                       "Cytokine;Interleukin"                           
#[17] "Cytokine;Rtk"                                    "Cytokine;Tgf"                                   
#[19] "Cytokine;Tnf family"                             "Glycoprotein"                                   
#[21] "Growth factor"                                   "Hormone"                                        
#[23] "Notch signalling"                                "Semaphorin"                                     
#[25] "Vasoconstrictor peptide"                         "Wnt family"                                     
#[27] "X"                                              
icellnet_subfamilies <- levels(as.factor(annotations$ICELLNET$subfamily))
#[1] "Adhesion molecule"       "Complement"              "Extracellular matrix"   
#[4] "Gla-containing protein"  "Hormone"                 "IL1."                   
#[7] "IL17"                    "RTK"                     "Semaphorin"             
#[10] "TGF"                     "TNF"                     "type 1"                 
#[13] "type 2"                  "Vasoconstrictor peptide" "Wnt family" 

italk_subclasses <- levels(as.factor(annotations$iTALK$subclass))
#[1] "checkpoint"    "cytokine"      "growth factor" "other" 

#####
# Location info:
#####

baccin2019_location <- levels(as.factor(annotations$Baccin2019$location))
#[1] "both"     "ecm"      "membrane" "other"    "secreted"
connectomedb2020_location <- levels(as.factor(annotations$connectomeDB2020$location))
#[1] "ECM"             "plasma membrane" "secreted"   

# Pathway info:
cellchatdb_pathways <- levels(as.factor(annotations$CellChatDB$pathway))

