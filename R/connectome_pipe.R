#' Function to call connectome with databases from OmniPath
#'
#' @param op_resource OmniPath Intercell Resource DN
#' @param sce Seurat object as input
#' @param .format bool whether to format output
#' @param ... dot params passed to connectome
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
#' @importFrom Seurat ScaleData
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr arrange select mutate distinct
#'
#' @export
call_connectome <- function(sce,
                            op_resource = NULL,
                            .format = TRUE,
                            ...){

    # Convert sce to seurat
    if(class(sce) == "SingleCellExperiment"){
        sce %<>% .liana_convert(., assay=assay)
    }

    if(!is.null(op_resource)){

        lr_db <- conn_formatDB(op_resource)
        lr_symbols <- union(lr_db$source_genesymbol,
                            lr_db$target_genesymbol)

        conn <- .conn_create(sce,
                             lr_symbols = lr_symbols,
                             lr_db,
                             ...
                             )

    } else{

        lr_db <- Connectome::ncomms8866_human %>%
            filter(Pair.Evidence == "literature supported")

        lr_symbols = union(lr_db$Ligand.ApprovedSymbol,
                           lr_db$Receptor.ApprovedSymbol)

        conn <- .conn_create(sce,
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
#' @import tibble
#'
#' @export
FormatConnectome <- function(conn,
                             ...){
    conn <- conn %>%
        Connectome::FilterConnectome(remove.na=TRUE) %>%
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
#' @param ... arguments passed `CreateConnectome` from `Connectome`
#'
#' @noRd
.conn_create <- function(sce,
                        lr_symbols,
                        lr_db,
                        ...){

    filt_genes <- lr_symbols[lr_symbols %in% rownames(sce)]
    sce <- Seurat::ScaleData(object = sce,
                                       features = filt_genes)

    conn <- Connectome::CreateConnectome(sce,
                                         LR.database = 'custom',
                                         custom.list = lr_db,
                                         ...
                                         )
    return(conn)
}
