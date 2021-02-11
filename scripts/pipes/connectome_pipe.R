#' Function to call connectome with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param .format bool whether to format output
#' @inheritDotParams Connectome::CreateConnectome
#' @return An unfiltered connectome df
#'
#' @details
#' Stats:
#' 1) The ‘weight_norm’ edge attribute is derived from the normalized expression
#'  of the ligand and the receptor in the single-cell data.
#' 2) The ‘weight_scale’ edge attribute is derived from the z-scores of the ligand
#'  and the receptor in each edge, and is of higher value when the ligand and receptor
#'   are more specific to a given pair of cell types
#' 3) p-val
call_connectome <- function(op_resource,
                            seurat_object,
                            .format = TRUE,
                            ...){
    library(Connectome)

    if(!is.null(op_resource)){

        # Format db to connectome
        lr_db <- op_resource %>%
            select("source_genesymbol", "target_genesymbol") %>%
            mutate(mode = "UNCAT") %>% # mode refers to interaction categories
            arrange(.$source_genesymbol) %>%
            as.data.frame()


        # scale genes to ligands and receptors available in the resource
        connectome.genes <- union(lr_db$source_genesymbol, lr_db$target_genesymbol)
        genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]
        seurat_object <- ScaleData(seurat_object, features = genes)

        # create connectome
        conn <- CreateConnectome(seurat_object,
                                 custom.list = lr_db,
                                 ...)

    } else{
        connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                                  Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
        genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]
        seurat_object <- ScaleData(seurat_object, features = genes)


        # HERE NOTE p.values
        conn <- CreateConnectome(seurat_object,
                                 species = 'human',
                                 ...)

    }

    if(.format){
        # here I use parameters as in their comparison when comp to CellPhoneDB
        conn <- conn %>% FormatConnectome(conn = .,
                                          min.pct = 0.1,
                                          min.z = 1,
                                          max.p = 0.05,
                                          remove.na = TRUE)
    }


    return(conn)
}

#' Helper function to filter and format connectome
#' @param conn connectome object
#' @inheritDotParams connectome::FilterConnectome
FormatConnectome <- function(conn,
                             ...){
    conn <- conn %>%
        FilterConnectome(.,
                         ...) %>%
        select(source, target,
               ligand, receptor,
               weight_norm,
               # weight_sc,
               p_val_adj.lig,
               p_val_adj.rec) %>%
        mutate(weight_norm = log(.$weight_norm + 1.0)
               # weight_sc = log(.$weight_sc + 1.0)
        )
}

