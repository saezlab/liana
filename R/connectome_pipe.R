#' Function to call connectome with databases from OmniPath
#'
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
                            op_resource = NULL,
                            .format = TRUE,
                            ...){

    if(!is.null(op_resource)){

        lr_db <- conn_formatDB(op_resource)
        lr_symbols <- union(lr_db$source_genesymbol,
                            lr_db$target_genesymbol)

        conn <- .conn_create(seurat_object,
                            lr_symbols = lr_symbols,
                            lr_db,
                            ...
        )

    } else{

        lr_db <- Connectome::ncomms8866_human %>%
            filter(Pair.Evidence == "literature supported")

        lr_symbols = union(lr_db$Ligand.ApprovedSymbol,
                           lr_db$Receptor.ApprovedSymbol)


        conn <- .conn_create(seurat_object,
                            lr_symbols = lr_symbols,
                            lr_db,
                            ...
                            )
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



#' Helper function to create a conn object
#'
#' @inheritParams call_connectome
#' @param lr_symbols ligand-receptor gene symbols
#' @param lr_db ligand-receptor resource
#' @inheritDotParams  Connectome::CreateConnectome
#'
#' @import Connectome
#'
#' @noRd
.conn_create <- function(seurat_object,
                        lr_symbols,
                        lr_db,
                        ...){
    filt_genes <- lr_symbols[lr_symbols %in% rownames(seurat_object)]
    seurat_object <- Seurat::ScaleData(object = seurat_object,
                                       features = filt_genes)

    conn <- CreateConnectome(seurat_object,
                             LR.database = 'custom',
                             custom.list = lr_db,
                             ...
                             )
    return(conn)
}
