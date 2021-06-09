#' Function to call SingleCellSignalR with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @param assay Seurat assay data to use
#'
#' @return An unfiltered iTALK df sorted by relevance
#' @importFrom Seurat GetAssayData Idents
#' @import SCAomni
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr distinct select
#'
#' @details
#' Stats:
#' LRScore = sqrt(LR product)/mean(raw counts) * sqrt(LR product) where
#' expression of l > 0 and r > 0
#' LRScore = 1 is the highest (~ most likely hit), 0 is the lowest.
#'
#' @export
call_sca <- function(op_resource,
                     seurat_object,
                     .format = TRUE,
                     assay = "SCT",
                     ...) {
  # Format OmnipathR resource
  if(!is.null(op_resource)){
    op_resource %<>% sca_formatDB
  } else{
    if(file.exists("input/LRdb.rda")){
      load("input/LRdb.rda")
      op_resource <- LRdb
    }
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
#' @importFrom purrr pluck
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#'
#' @export
FormatSCA <- function(sca_res, remove.na = TRUE) {
  sca_res <- sca_res %>%
    pluck("full-network") %>%
    separate(ligand,
             into = c("source", "ligand"),
             sep = "[.]") %>%
    separate(receptor,
             into = c("target", "receptor"),
             sep = "[.]") %>%
    select(source, ligand, target, receptor, LRscore) %>%
    as_tibble()
  return(sca_res)
}


#' Helper Function to convert Omni to LRdb Format
#' @param op_resource OmniPath resource
#' @export
sca_formatDB <- function(op_resource){
  op_resource %>%
  select(ligand = source_genesymbol,
         receptor = target_genesymbol,
         source = sources,
         PMIDs = references) %>%
    distinct()
}
