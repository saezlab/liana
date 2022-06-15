#' Run iTALK with OmniPath data [[DEPRECATED]]
#'
#' @param sce Seurat object or SingleCellExperiment as input
#' @param op_resource OmniPath Intercell Resource DN
#' @param assay assay to use from Seurat object
#' @param .format bool: whether to format output
#' @param ... Parameters passed to Seurat FindMarkers (ref requires import)
#'
#' @return An unfiltered iTALK df sorted by relevance
#'
#' In this case, we use the product of the logFC rather than thresholding, as in
#' the original implementation.
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unite expand_grid
#' @importFrom dplyr select rename mutate group_by group_split
#'
#' @export
#'
#' @details In order to be comparable with the remainder of the methods, we
#' calculate the mean of the ligand and receptor logFC.
#' The original implementation only uses the DE genes above a certain logFC
#' threshold.
call_italk <- function(
  sce,
  op_resource,
  assay = "RNA",
  .format = TRUE,
  ...){

  # Convert sce to seurat
  if(class(sce) == "SingleCellExperiment"){
    sce %<>% .liana_convert(., assay=assay)
  }

  if(!is.null(op_resource)){
    op_resource %<>% italk_formatDB
  }

  # create a dataframe of the cell labels
  cell_type <- Idents(sce) %>%
    data.frame(group = ., row.names = names(.)) %>%
    unname() %>%
    as.vector()


  deg <- Seurat::FindAllMarkers(sce,
                                assay = assay,
                                ...
                                ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split()  %>%
    map(function(x){

      logcol <- ifelse("avg_logFC" %in% colnames(x), # Seurat...
                       "avg_logFC", # version 3.2.3
                       "avg_log2FC") # version 4.0.3

      x %>%
        rename(p.value = 'p_val',
               logFC = !!logcol,
               q.value = 'p_val_adj',
               cell_type = 'cluster',
               gene = 'gene')
      }) %>%
    setNames(levels(Idents(sce)))

  # Iterate over cell type pairs and Find LR
  idents <- as.character(unique(Idents(sce)))
  comb <- expand_grid(source = idents, target = idents)

  res <- comb %>%
    pmap(function(source, target){
      iTALK::FindLR(
        deg[[source]],
        deg[[target]],
        datatype = 'DEG',
        comm_type = 'other',
        database = op_resource
        )
      }) %>%
    bind_rows()

  if (.format) {
    res <- res %>% FormatiTALK(remove.na = TRUE)
  }

  return(res)
}




#' Helper Function to Filter and format iTalk results
#'
#' @param italk_res iTalk results object
#' @param remove.na bool whether to filter NA
#'
#' @importFrom tibble tibble
#' @export
FormatiTALK <- function(italk_res,
                        remove.na = TRUE){

    italk_res <- tibble(
      'source' = italk_res$cell_from,
      'ligand' = italk_res$ligand,
      'target' = italk_res$cell_to,
      'receptor' = italk_res$receptor,
      'logFC_from' = italk_res$cell_from_logFC,
      'logFC_to' = italk_res$cell_to_logFC,
      'qval_from' = italk_res$cell_from_q.value,
      'qval_to' = italk_res$cell_to_q.value
      ) %>%
      mutate(logfc_comb = mean(c(logFC_from, logFC_to)))

  return(italk_res)
}


#' Helper Function to convert Omni to iTALK resource Format
#'
#' @param op_resource OmniPath resource
#' @export
#'
#' @noRd
italk_formatDB <- function(op_resource){
  op_resource %>%
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
