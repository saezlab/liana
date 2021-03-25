#' Function to call iTalk with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param assay assay to use from Seurat object
#' @param .format bool: whether to format output
#' @param .DE bool: whether to use DE (TRUE) or highlyVarGenes (FALSE)
#' @param .deg if is NULL run FindAllMarkers
#' @inheritDotParams iTALK::rawParse
#' @return An unfiltered iTALK df sorted by relevance
#'
#' @details
#' Stats:
#' Ligand and Receptor Expressions (and P-values if ran with .DE==TRUE)
#' Dot params are inherited from Seurat::FindAllMarkers, if .deg = TRUE
#' @import iTALK Seurat tidyverse
#' @export
call_italk <- function(op_resource,
                       seurat_object,
                       assay = "SCT",
                       .format = TRUE,
                       .DE = FALSE,
                       .deg = NULL,
                       ...){

  if(!is.null(op_resource)){
    op_resource <- op_resource %>%
      unite(col = "Pair.Name", source_genesymbol, target_genesymbol,
            sep="_", remove = FALSE) %>%
      rename('Ligand.ApprovedSymbol' = source_genesymbol,
             'Receptor.ApprovedSymbol' = target_genesymbol) %>%
      mutate("Classification" = "other",
             'Receptor.Name' = Receptor.ApprovedSymbol,
             'Ligand.Name' = Ligand.ApprovedSymbol) %>%
      select(Pair.Name, Ligand.ApprovedSymbol, Ligand.Name,
             Receptor.ApprovedSymbol, Receptor.Name, Classification) %>%
      as.data.frame()
  }

  # create a dataframe of the cell labels
  cell_type <- Idents(seurat_object) %>%
    data.frame(group = ., row.names = names(.)) %>%
    unname() %>%
    as.vector()

  # Prepare data from seurat object
  input_data <-
    GetAssayData(seurat_object, assay = assay, slot = "data") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    cbind(., cell_type) # bind cell types from meta


  if(.DE){
    if (is.null(.deg)) {
      deg <- FindAllMarkers(seurat_object,
                            assay = assay,
                            ...) %>%
        dplyr::group_by(cluster) %>%
        dplyr::group_split() %>%
        map(function(x)
          x %>%
            rename(p.value = 'p_val',
                   logFC = 'avg_logFC',
                   pct.1 = 'pct.1',
                   pct.2 = 'pct.2',
                   q.value = 'p_val_adj',
                   cell_type = 'cluster',
                   gene = 'gene')) %>%
        setNames(levels(Idents(seurat_object)))
    }

    # Iterate over cell type pairs and Find LR
    comb <- combn(levels(Idents(seurat_object)), 2)
    res = list()
    for (i in  seq_len(dim(comb)[2])) {
      res[[paste0(comb[, i][1], '_x_', comb[, i][2])]] <-
        FindLR(deg[[comb[, i][1]]], deg[[comb[, i][2]]],
               datatype = 'DEG',
               comm_type = 'other',
               database = op_resource)
    }
    res <- bind_rows(res)

  } else{
    highly_exprs_genes <- rawParse(input_data,
                                   ...)
    res <- FindLR(
      highly_exprs_genes,
      datatype = 'mean count',
      comm_type = 'other',
      database = op_resource
    ) %>%
      arrange(desc(cell_from_mean_exprs * cell_to_mean_exprs))
  }

  if (.format) {
    res <- res %>% FormatiTALK(remove.na = TRUE, .DE = .DE)
  }

  return(res)
}




#' Helper function to filter and format iTalk results
#' @param italk_res iTalk results object
#' @param remove.na bool whether to filter NA
#' @export
FormatiTALK <- function(italk_res,
                        remove.na = TRUE,
                        .DE = FALSE){
  if(!.DE){
    italk_res <- tibble(
      'source' = italk_res$cell_from,
      'ligand' = italk_res$ligand,
      'target' = italk_res$cell_to,
      'receptor' = italk_res$receptor,
      'weight_from' = italk_res$cell_from_mean_exprs,
      'weight_to' = italk_res$cell_to_mean_exprs
    ) %>%
      mutate(weight_comb = abs(weight_from * weight_to))
    if (remove.na) {
      italk_res <- italk_res[!(is.na(italk_res$weight_from) &
                                 is.na(italk_res$weight_to)), ]
    }

  } else{
    italk_res <- tibble(
      'source' = italk_res$cell_from,
      'ligand' = italk_res$ligand,
      'target' = italk_res$cell_to,
      'receptor' = italk_res$receptor,
      'logFC_from' = italk_res$cell_from_logFC,
      'logFC_to' = italk_res$cell_to_logFC,
      'qval_from' = italk_res$cell_from_q.value,
      'qval_to' = italk_res$cell_to_q.value,
    )  %>%
      mutate(weight_comb = abs(logFC_from * logFC_to))
  }

  return(italk_res)
}
