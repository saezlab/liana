#' Function to call SingleCellSignalR with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @param assay
#' @param genes_cutoff highly variable genes top_n cutoff
#' @param stats stats considered to select top_n
#' @return An unfiltered iTALK df sorted by relevance
#'
#' @details
#' Stats:
#' 1) The ‘weight_norm’ edge attribute is derived from the normalized expression
#'  of the ligand and the receptor in the single-cell data.
#' 2) The ‘weight_scale’ edge attribute is derived from the z-scores of the ligand
#'  and the receptor in each edge, and is of higher value when the ligand and receptor
#'   are more specific to a given pair of cell types
#' 3) p-val
call_sca <- function(op_resource,
                       seurat_object,
                      .format = TRUE,
                       assay = 'SCT',
                            ...){
    require(Seurat)
    require(SCAomni)
    require(dplyr)
    # Selecting the OmnipathR resource
    sel <- op_resource[,c('source_genesymbol','target_genesymbol','sources','references')]
    sel <- sel %>% rename(c(ligand='source_genesymbol',receptor='target_genesymbol',source='sources',PMIDs='references'))
    # Preparing data from seurat object
    # normalized data matrix
    data.input <- GetAssayData(seurat_object, assay = assay, slot = "data")
    labels <- Idents(seurat_object)
    # create a dataframe of the cell labels
    meta <- data.frame(group = labels, row.names = names(labels))
    data <- t(as.matrix(data.input))
    signal <- cell_signaling(data=data1,genes=row.names(count_data),cluster=as.numeric(meta_data),c.names = levels(meta_data),s.score = 0.05,species='homo sapiens',LRdb = sel,write = F)
    res <- inter_network(data = data1, signal = signal,genes=row.names(count_data),cluster=as.numeric(meta_data),c.names = levels(meta_data), write = FALSE)
    if(.format){
        res <- res %>% FormatSCA(.data)
    }
    return(res)
}


#' Function to call SingleCellSignalR with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @param assay
#' @param .deg if is NULL run FindAllMarkers
#' @param stats stats considered to select top_n
#' @return An unfiltered iTALK df sorted by relevance
#'
#' @details
#' Stats:
#' 1) The ‘weight_norm’ edge attribute is derived from the normalized expression
#'  of the ligand and the receptor in the single-cell data.
#' 2) The ‘weight_scale’ edge attribute is derived from the z-scores of the ligand
#'  and the receptor in each edge, and is of higher value when the ligand and receptor
#'   are more specific to a given pair of cell types
#' 3) p-val
call_sca_de <- function(op_resource,
                       seurat_object,
                      .format = TRUE,
                       assay = 'SCT',
                       .deg = NULL,
                            ...){
    require(Seurat)
    require(SCAomni)
    require(dplyr)
    # Selecting the OmnipathR resource
    sel <- op_resource[,c('source_genesymbol','target_genesymbol','sources','references')]
    sel <- sel %>% rename(c(ligand='source_genesymbol',receptor='target_genesymbol',source='sources',PMIDs='references'))
    # Preparing data from seurat object
    # normalized data matrix
    data.input <- GetAssayData(seurat_object, assay = assay, slot = "data")
    labels <- Idents(seurat_object)
    # create a dataframe of the cell labels
    meta <- data.frame(group = labels, row.names = names(labels))
    data <- t(as.matrix(data.input))
    signal <- cell_signaling(data=data1,genes=row.names(count_data),cluster=as.numeric(meta_data),c.names = levels(meta_data),s.score = 0.05,species='homo sapiens',LRdb = sel,write = F)
    res <- inter_network(data = data1, signal = signal,genes=row.names(count_data),cluster=as.numeric(meta_data),c.names = levels(meta_data), write = FALSE)
    return(res)
}

#' Helper function to filter and format connectome
#' @param result object
#' @param remove.na
FormatSCA <- function(conn,remove.na=TRUE){
      conn <- conn$`full-network` %>% separate(ligand,c("Ligand.Cell","Ligand"),sep="[.]")
      conn <- conn %>% separate(receptor,c("Receptor.Cell","Receptor"),sep="[.]")
      names(conn)[length(names(sca_res))] <- 'weight'
      return(conn)
}
