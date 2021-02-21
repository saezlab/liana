#' Function to call SingleCellSignalR with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @param assay
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
                     assay = "SCT",
                     .default_db = FALSE,
                     ...) {
  require(Seurat)
  require(SCAomni)
  require(dplyr)

  # Format OmnipathR resource
  if(!.default_db){
    op_resource <- op_resource %>%
      select(ligand = source_genesymbol,
             receptor = target_genesymbol,
             source = sources,
             PMIDs = references)
  }

  # Prepare data from Seurat object
  input_data <-
    GetAssayData(seurat_object, assay = assay, slot = "data") %>%
    as.matrix()
  labels <- Idents(seurat_object)

  # Compute interactions between cell clusters
  signal <- cell_signaling(data = input_data,
                           genes = row.names(input_data),
                           cluster = as.numeric(labels),
                           c.names = levels(Idents(seurat_object)),
                           species = 'homo sapiens',
                           LRdb = op_resource,
                           write = FALSE,
                           ...
    )

  # Compute intercellular gene networks
  sca_res <- inter_network(data = input_data,
                           signal = signal,
                           genes = row.names(input_data),
                           cluster = as.numeric(labels),
                           c.names = levels(Idents(seurat_object)),
                           write = FALSE
    )


  if (.format) {
    sca_res <- sca_res %>% FormatSCA(.data)
  }
  return(sca_res)
}



#' Helper function to format SingleCellSignalR results
#' @param sca_res Unformatted SCA results
#' @param remove.na bool whether to filter SCA output
FormatSCA <- function(sca_res, remove.na = TRUE) {
  sca_res <- sca_res %>%
    pluck("full-network") %>%
    separate(ligand,
             into = c("source", "ligand"),
             sep = "[.]") %>%
    separate(receptor,
             into = c("target", "receptor"),
             sep = "[.]") %>%
    select(source, ligand, target, receptor, LRscore)
  return(sca_res)
}
