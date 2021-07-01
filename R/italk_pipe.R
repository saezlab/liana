#' Run iTALK with OmniPath data
#'
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param assay assay to use from Seurat object
#' @param .format bool: whether to format output
#' @param .DE bool: whether to use DE (TRUE) or highlyVarGenes (FALSE)
#' @inheritDotParams Seurat::FindAllMarkers
#'
#' @return An unfiltered iTALK df sorted by relevance
#'
#' @details
#' Stats:
#' Ligand and Receptor Expressions (and P-values if ran with .DE==TRUE)
#' Dot params are inherited from Seurat::FindAllMarkers, if .deg = TRUE
#' @importFrom iTALK FindLR
#' @importFrom Seurat Idents FindAllMarkers GetAssayData
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr unite expand_grid
#' @importFrom dplyr select rename mutate group_by group_split
#'
#' @export
call_italk <- function(
    op_resource,
    seurat_object,
    assay = "RNA",
    .format = TRUE,
    .DE = TRUE,
    ...
){

  if(!is.null(op_resource)){
    op_resource %<>% italk_formatDB
  }

  # create a dataframe of the cell labels
  cell_type <- Idents(seurat_object) %>%
    data.frame(group = ., row.names = names(.)) %>%
    unname() %>%
    as.vector()


  deg <- Seurat::FindAllMarkers(seurat_object,
                                assay = assay,
                                ...
                                ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split() %>%
    map(function(x){
      x %>%
        rename(p.value = 'p_val',
               logFC = 'avg_logFC',
               q.value = 'p_val_adj',
               cell_type = 'cluster',
               gene = 'gene')
      }) %>%
    setNames(levels(Idents(seurat_object)))

  # Iterate over cell type pairs and Find LR
  idents <- as.character(unique(Idents(seurat_object)))
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
    res <- res %>% FormatiTALK(remove.na = TRUE, .DE = .DE)
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
                        remove.na = TRUE,
                        .DE = FALSE){

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
      mutate(weight_comb = (logFC_from * logFC_to))

  return(italk_res)
}


#' Helper Function to convert Omni to iTALK resource Format
#' @param op_resource OmniPath resource
#' @export
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
