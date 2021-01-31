# Dependancies

# library(devtools)
# devtools::install_github("https://github.com/renozao/NMF")
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("sqjin/CellChat")
library(CellChat)
library(tidyverse)
library(patchwork)
library(ggalluvial)
library(igraph)
library(tidyverse)
library(OmnipathR)
options(stringsAsFactors = FALSE)

# library(Seurat)
# library(OmnipathR)

seurat_object <- readRDS("output/pbmc3k_processed.rds")


# Create cellchat object
cellchat <- createCellChat(object = seurat_object,
                           group.by = "ident")



# data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
# labels <- Idents(seurat_object)
# meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
#
# unique(meta$labels)

# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


## use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

## set the used database in the object
cellchat@DB <- CellChatDB.use


## subset the expression data of signaling genes
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) # do parallel

# To infer the cell state-specific communications,
# we identify over-expressed ligands or receptors in one cell group and then
# identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## Compute the communication probability and infer cellular communication network
# i.e. project gene expression data onto protein-protein interaction (PPI) network.
# Specifically, a diffusion process is used to smooth genes’ expression values
# based on their neighbors’ defined in a high-confidence
# experimentally validated protein-protein network
# cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)
head(df.net)


## Infer the cell-cell communication at a signaling pathway level
# CellChat computes the communication probability on signaling pathway level
# by summarizing the communication probabilities of all ligands-receptors
# interactions associated with each signaling pathway.

cellchat <- computeCommunProbPathway(cellchat)
# NB: The inferred intercellular communication network of each ligand-receptor
# pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.






# Connectome_Omni Pipe ----
seurat_object <- readRDS("output/pbmc3k_processed.rds")

cellchat.omni <- createCellChat(object = seurat_object,
                                group.by = "ident")

# load CellChatDB
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data


# OmniPath
op_resource <- omni_resources[["talklr"]]

# Overwrite with Omni
interaction_input <- CellChatDB$interaction

# Compare OP vs CellChat format
head(op_resource)
head(interaction_input)


# get complexes and interactions from omnipath
complex_interactions <- op_resource %>%
    select("ligand" = source_genesymbol,
           "receptor" = target_genesymbol,
           "evidence" = sources,
           category_intercell_source,
           category_intercell_target,
           genesymbol_intercell_source,
           genesymbol_intercell_target,
           is_directed,
           is_stimulation,
           is_inhibition
           ) %>%
    unite("annotation",
          c(category_intercell_source, category_intercell_target),
          sep="-") %>%
    # join with CellChat resource for additional fields
    left_join(interaction_input %>%
                  select(-c(evidence, annotation))) %>%
    unite("interaction_name", c(ligand, receptor), remove = FALSE) %>%
    mutate_at(vars(everything()), ~ replace(., is.na(.), ""))


# Get OmniPath directed info
omni_directions <- complex_interactions %>%
    select(interaction_name,
           is_directed,
           is_stimulation,
           is_inhibition,
           co_A_receptor,
           co_I_receptor) %>%
    mutate(direction = pmap(., function(interaction_name,
                                is_directed,
                                is_stimulation,
                                is_inhibition,
                                co_A_receptor,
                                co_I_receptor){
        suppressMessages(message(interaction_name))
        if(co_A_receptor=="" & co_I_receptor==""){
            if(is_directed==1){
                if(is_stimulation==1 & is_inhibition==1){
                    return("Both")
                } else if(is_stimulation==1){
                    return("Stimulation")
                } else if(is_inhibition==1){
                    return("Inhibition")
                }
            }
        } else{
            return(NA)
        }
    })) %>% unnest(direction)


omni_interactions <- complex_interactions %>%
    select(all_of(names(interaction_input))) %>%
    # set LR interactions as rowname
    left_join(., (omni_directions %>%
                      select(interaction_name,
                             direction))) %>%
    mutate_at(vars(everything()), ~ replace(., is.na(.), "")) %>%
    mutate(co_A_receptor = ifelse(.data$co_A_receptor == "" & (direction == "Stimulation") | (direction == "both"),
                                  "Stimulation", .data$co_A_receptor),
           co_I_receptor = ifelse(.data$co_I_receptor == "" & (direction == "Inhibition") | (direction == "both"),
                                  "Inhibition", .data$co_I_receptor)) %>%
    mutate("interaction_name2" = interaction_name) %>%
        distinct_at(.vars="interaction_name2", .keep_all = TRUE) %>%
        column_to_rownames("interaction_name2")



# show summaries
summary(as.factor(interaction_input$annotation))
summary(as.factor(omni_interactions$annotation))




# Get Omni Complexes
omni_complexes <- complex_interactions %>%
    filter(str_detect(genesymbol_intercell_source, "COMPLEX") |
               str_detect(genesymbol_intercell_target, "COMPLEX")) %>%
    select(genesymbol_intercell_source,
           genesymbol_intercell_target)

# Convert to CellChat format
omni_complexes <- union(omni_complexes$genesymbol_intercell_source,
                        omni_complexes$genesymbol_intercell_target) %>%
    str_subset(pattern = "COMPLEX") %>%
    str_replace(pattern = "COMPLEX:", "") %>%
    enframe() %>%
    separate(col=value, sep="_",
             into = c("subunit_1", "subunit_2",
                      "subunit_3", "subunit_4"), remove=FALSE) %>%
    mutate_at(vars(everything()), ~ replace(., is.na(.), "")) %>%
    select(-name) %>%
    column_to_rownames("value")





# Replace Default DB with OmniPath Resource
# Here, I filter ECM, as when I don't an error is encountered
# both with mine and their dataset (they too filter for secreted signaling only)
CellChatDB.omni <- CellChatDB
CellChatDB.omni$interaction <- complex_interactions %>%
    select(-c(genesymbol_intercell_source,
              genesymbol_intercell_target)) %>%
    filter(annotation!="ecm-receptor")

CellChatDB.omni$complex <- omni_complexes


## set the used database in the object
cellchat.omni@DB <- CellChatDB.omni

## use a subset of CellChatDB for cell-cell communication analysis

## set the used database in the object

## subset the expression data of signaling genes
cellchat.omni <- subsetData(cellchat.omni)
future::plan("multiprocess", workers = 4) # do parallel

# To infer the cell state-specific communications,
# we identify over-expressed ligands or receptors in one cell group and then
# identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
cellchat.omni <- identifyOverExpressedGenes(cellchat.omni)
cellchat.omni <- identifyOverExpressedInteractions(cellchat.omni)



## Compute the communication probability and infer cellular communication network
# i.e. project gene expression data onto protein-protein interaction (PPI) network.
# Specifically, a diffusion process is used to smooth genes’ expression values
# based on their neighbors’ defined in a high-confidence
# experimentally validated protein-protein network
cellchat.omni <- projectData(cellchat.omni, PPI.human)
cellchat.omni <- computeCommunProb(cellchat.omni,
                                   raw.use=T)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.omni <- filterCommunication(cellchat.omni,
                                     min.cells = 10)


# Extract the inferred cellular communication network as a data frame
df.omni <- subsetCommunication(cellchat.omni,
                               thresh= 0.05)
?subsetCommunication

head(df.omni)

# !!! It is worth noting that our results seem to overestimate the LR hits
# from the CellChat Package - this is likely due to differences into the
# curation effort that went into identifying co-inhibitors and co-activators
# Information about how and which fields of the DB are used is lacking
# and as such without further information it would be difficult to reproduce
# I submitted an issue to their GH to provide more info about the DB

# Maybe we keep its won dataset and benchmark it only with it versus the other tools
