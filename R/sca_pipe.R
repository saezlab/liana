#' Function to call SingleCellSignalR with databases from OmniPath [[DEPRECATED]]
#'
#' @param sce SingleCellExperiment or SeuratObject as input
#' @param op_resource OmniPath Intercell Resource DN
#' @param .format bool whether to format output
#' @param assay Seurat assay data to use
#' @param assay.type count slot (logcounts by default)
#' @param ... arguments passed to `SCAomni::cell_signaling`

#' @importFrom SeuratObject GetAssayData Idents
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
#'
#' @return An unfiltered SCA tibble
call_sca <- function(sce,
                     op_resource,
                     .format = TRUE,
                     assay = "RNA",
                     assay.type = "logcounts",
                     ...){

  # Convert sce to seurat
  if(class(sce) == "SingleCellExperiment"){
    sce %<>% .liana_convert(., assay=assay)
  }

  if(class(sce) == "Seurat" & assay.type=="logcounts"){
    assay.type = "data"
  }

  # Format OmnipathR resource
  if(!is.null(op_resource)){
    op_resource %<>% sca_formatDB
  } else{
    if(file.exists(system.file(package = "liana", "LRdb.rda"))){
      load(system.file(package = "liana", "LRdb.rda"))
      op_resource <- LRdb
    } else{
      stop("Could not locate LRdb.rda")
    }
  }

  # Prepare data from Seurat object
  input_data <-
    SeuratObject::GetAssayData(sce,
                         assay = assay,
                         slot = assay.type)
  labels <- SeuratObject::Idents(sce)

  # Compute interactions between cell clusters
  signal <- SCAomni::cell_signaling(data = input_data,
                                    genes = row.names(input_data),
                                    cluster = as.numeric(labels),
                                    c.names = levels(Idents(sce)),
                                    species = 'homo sapiens',
                                    LRdb = op_resource,
                                    int.type="autocrine", # includes both para and auto...
                                    write = FALSE,
                                    verbose = FALSE,
                                    ...
                                    )

  # Compute intercellular gene networks
  sca_res <- SCAomni::inter_network(data = input_data,
                                    signal = signal,
                                    genes = row.names(input_data),
                                    cluster = as.numeric(labels),
                                    c.names = levels(Idents(sce)),
                                    write = FALSE
                                    )
  if (.format) {
    sca_res %<>% FormatSCA
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
             sep = "⊎") %>%
    separate(receptor,
             into = c("target", "receptor"),
             sep = "⊎") %>%
    select(source, ligand, target, receptor, LRscore) %>%
    as_tibble()
  return(sca_res)
}


#' Helper Function to convert Omni to LRdb Format
#'
#' @param op_resource OmniPath resource
#'
#' @export
sca_formatDB <- function(op_resource){
  op_resource %>%
  select(ligand = source_genesymbol,
         receptor = target_genesymbol,
         source = sources,
         PMIDs = references) %>%
    distinct()
}
