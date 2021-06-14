#' Function to call connectome with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @inheritDotParams Connectome::CreateConnectome
#'
#' @return An unfiltered connectome results df
#'
#' @details
#' Stats:
#' 1) The ‘weight_norm’ edge attribute is derived from the normalized expression
#'  of the ligand and the receptor in the single-cell data.
#' 2) The ‘weight_scale’ edge attribute is derived from the z-scores of the ligand
#'  and the receptor in each edge, and is of higher value when the ligand and receptor
#'   are more specific to a given pair of cell types
#' 3) DEG p-values for L and R
#'
#' @import Connectome
#' @importFrom Seurat ScaleData
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr arrange select mutate distinct
#'
#' @export
call_connectome <- function(seurat_object,
                            op_resource,
                            .format = TRUE,
                            .spatial = FALSE,
                            ...){

    if(.spatial){
        seurat_object@assays$RNA <- seurat_object@assays$Spatial
        seurat_object@assays$Spatial <- NULL
    }

    if(!is.null(op_resource)){
        # Format db to connectome
        lr_db <- conn_formatDB(op_resource)

        # scale genes to ligands and receptors available in the resource
        connectome.genes <- union(lr_db$source_genesymbol, lr_db$target_genesymbol)
        genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]
        seurat_object <- ScaleData(seurat_object, features = genes)

        # create connectome
        conn <- CreateConnectome(seurat_object,
                                 LR.database = 'custom',
                                 custom.list = lr_db,
                                 ...)

    } else{
        lr_db <- Connectome::ncomms8866_human %>%
            filter(Pair.Evidence == "literature supported")
        connectome.genes <- union(lr_db$Ligand.ApprovedSymbol,
                                  lr_db$Receptor.ApprovedSymbol)
        genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]

        seurat_object <- ScaleData(object = seurat_object,
                                   features = genes)
        conn <- CreateConnectome(seurat_object,
                                 LR.database = 'custom',
                                 custom.list = lr_db,
                                 ...)
    }

    if(.format){
        conn %<>% FormatConnectome
    }

    return(conn)
}


#' Helper function to filter and format connectome
#'
#' @param conn connectome object
#' @importFrom Connectome FilterConnectome
#' @import tibble
#'
#' @export
FormatConnectome <- function(conn,
                             ...){
    conn <- conn %>%
        FilterConnectome(remove.na=TRUE) %>%
        select(source, target,
               ligand, receptor,
               weight_norm,
               weight_sc,
               p_val_adj.lig,
               p_val_adj.rec) %>%
        as_tibble()

}


#' Helper Function to convert Omni to Connectome resource Format
#' @param op_resource OmniPath resource
#' @export
conn_formatDB <- function(op_resource){
    op_resource %>%
        select("source_genesymbol", "target_genesymbol") %>%
        mutate(mode = "UNCAT") %>% # mode refers to interaction categories
        arrange(.$source_genesymbol) %>%
        distinct() %>%
        as.data.frame()
}
