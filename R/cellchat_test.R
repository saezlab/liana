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



#  CellChat x OmniPath

# Connectome_Omni Pipe ----
# OmniPath
op_resource <- omni_resources[["connectomeDB2020"]]

# Overwrite with Omni
interaction_input <- CellChatDB$interaction

# Compare OP vs CellChat format
head(op_resource)
head(interaction_input)

cellchat_cols <- names(interaction_input)

xx_interactions <-  op_resource %>%
    select("ligand" = source_genesymbol,
           "receptor" = target_genesymbol,
                  "evidence" = sources,
                  category_intercell_source,
                  category_intercell_target
                  ) %>%
           unite("annotation",
                 c(category_intercell_source, category_intercell_target),
                 sep="-") %>%
    left_join(interaction_input %>%
                  select(-c(evidence, annotation))) %>%
    unite("interaction_name", c(ligand, receptor), remove = FALSE) %>%
    mutate_at(vars(everything()), ~ replace(., is.na(.), "")) %>%
    select(all_of(cellchat_cols)) %>%
    as.data.frame()



# Compare now
head(xx_interactions)
head(interaction_input)


# Run with Omni DB
CellChatDB.omni <- CellChatDB
CellChatDB.omni$interaction <- xx_interactions


summary(as.factor(interaction_input$annotation))
summary(as.factor(CellChatDB.omni$interaction$annotation))



# Create cellchat object
cellchat.omni <- createCellChat(object = seurat_object,
                           group.by = "ident")




## use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use.omni <- subsetDB(CellChatDB.omni, "ligand-receptor")
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

## set the used database in the object
cellchat.omni@DB <- CellChatDB.use.omni


## subset the expression data of signaling genes
cellchat.omni <- subsetData(cellchat.omni)
# future::plan("multiprocess", workers = 4) # do parallel

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
# cellchat.omni <- projectData(cellchat.omni, PPI.human)
cellchat.omni <- computeCommunProb(cellchat.omni, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.omni <- filterCommunication(cellchat.omni, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.omni <- subsetCommunication(cellchat.omni)

head(df.omni)










# Update
# CellChatDB <- list()
# CellChatDB$interaction <- interaction_input
# CellChatDB$complex <- complex_input
# CellChatDB$cofactor <- cofactor_input
# CellChatDB$geneInfo <- geneInfo


