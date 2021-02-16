#' Function to call iTalk with databases from OmniPath
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
call_italk <- function(op_resource,
                       seurat_object,
                      .format = TRUE,
                       assay = 'SCT',
                       genes_cutoff = 50,
                       stats =  'mean',
                            ...){
    require(iTALK)
    # Selecting the OmnipathR resource
    sel <- op_resource
    dataset<-data.frame(
                 'Pair.Name' = paste0(sel$source_genesymbol,'_',sel$target_genesymbol),
                 'Ligand.ApprovedSymbol'=sel$source_genesymbol,
                 'Ligand.Name'=sel$source_genesymbol,
                 'Receptor.ApprovedSymbol'=sel$target_genesymbol,
                 'Receptor.Name'=sel$target_genesymbol,
                 'Classification'=rep('other',length(sel$target_genesymbol))
                )
    # Preparing data from seurat object
    # normalized data matrix
    data.input <- GetAssayData(seurat_object, assay = assay, slot = "data")
    labels <- Idents(seurat_object)
    # create a dataframe of the cell labels
    meta <- data.frame(group = labels, row.names = names(labels))
    data <- t(as.matrix(data.input))
    rownames(data) <- NULL
    names(meta) <- NULL
    cell_type = as.vector(meta)
    data <- as.data.frame(cbind(data,cell_type))
    highly_exprs_genes<-rawParse(data,
                                 top_genes=genes_cutoff,
                                 stats='mean')
    comm_list<-'other'
    res<-FindLR(highly_exprs_genes,
                datatype='mean count',
                comm_type=comm_list,
                database = as.data.frame(dataset))
    res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),]
    if(.format){
        # here I use parameters as in their comparison when comp to CellPhoneDB
        res <- res %>% FormatiTALK(conn = .,remove.na = TRUE)
    }
    return(res)
}


#' Function to call iTalk with databases from OmniPath
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
call_italk_de <- function(op_resource,
                       seurat_object,
                      .format = TRUE,
                       assay = 'SCT',
                       stats =  'mean',
                       .deg = NULL,
                            ...){
    require(iTALK)
    require(dplyr)
    # Selecting the OmnipathR resource
    sel <- op_resource
    dataset<-data.frame(
                 'Pair.Name' = paste0(sel$source_genesymbol,'_',sel$target_genesymbol),
                 'Ligand.ApprovedSymbol'=sel$source_genesymbol,
                 'Ligand.Name'=sel$source_genesymbol,
                 'Receptor.ApprovedSymbol'=sel$target_genesymbol,
                 'Receptor.Name'=sel$target_genesymbol,
                 'Classification'=rep('other',length(sel$target_genesymbol))
                )
    # Preparing data from seurat object
    # normalized data matrix
    data.input <- GetAssayData(seurat_object, assay = assay, slot = "data")
    labels <- Idents(seurat_object)
    # create a dataframe of the cell labels
    meta <- data.frame(group = labels, row.names = names(labels))
    data <- t(as.matrix(data.input))
    rownames(data) <- NULL
    names(meta) <- NULL
    cell_type = as.vector(meta)
    data <- as.data.frame(cbind(data,cell_type))
    if(is.null(.deg)){
        deg <- FindAllMarkers(seurat_object,assay ='RNA')
        deg <-deg %>%
              dplyr::group_by(cluster) %>%
              dplyr::group_split()
        deg <-lapply(deg, function(x)  rename(x, c(p.value='p_val',logFC='avg_logFC',pct.1='pct.1',pct.2='pct.2',q.value='p_val_adj',cell_type='cluster',gene='gene')))
        names(deg) <- levels(Idents(seurat_object))
    }
    comb <- combn(levels(Idents(seurat_object)),2)
    res_all = list()
    for( i in  seq_len(dim(comb)[2])){
      print(comb[,i])
      res_all[[paste0(comb[,i][1],'_x_',comb[,i][2])]]<-FindLR(deg[[comb[,i][1]]],deg[[comb[,i][2]]],datatype='DEG',comm_type='other')
    }
    res_all <- bind_rows(res_all)
    if(.format){
        # here I use parameters as in their comparison when comp to CellPhoneDB
        res_all <- res_all %>% FormatiTALK(conn = .,remove.na = TRUE)
    }
    return(res_all)
}

#' Helper function to filter and format connectome
#' @param conn connectome object
#' @param remove.na connectome object
FormatiTALK <- function(conn,remove.na=TRUE){
   conn <- tibble('Ligand.Cell'=conn$cell_from,
                 'Ligand'=conn$ligand,
                 'Receptor.Cell'=conn$cell_to,
                 'Receptor'=conn$receptor,
                 'weight_from'=conn$cell_from_mean_exprs,
                 'weight_to'=conn$cell_to_mean_exprs
                  )
    if(remove.na){
            conn <- conn[!(is.na(conn$weight_from) & is.na(conn$weight_to)),]
    }
}
