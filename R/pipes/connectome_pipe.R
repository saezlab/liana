#' Function to call connectome with databases from OmniPath
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @inheritDotParams Connectome::CreateConnectome
#' @return An unfiltered connectome df
call_connectome <- function(op_resource,
                            seurat_object,
                            ...){
    library(Connectome)

    op_resource <- omni_resources$connectomeDB2020

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
    return(conn)
}


#' Function to call connectome with default DB
#' @param seurat_object Seurat object as input
#' @inheritDotParams Connectome::CreateConnectome
#' @return An unfiltered Connectome df
call_connectome_default <- function(seurat_object,
                                    ...){
    connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                              Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
    genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]
    seurat_object <- ScaleData(seurat_object, features = genes)


    # HERE NOTE p.values
    conn <- CreateConnectome(seurat_object,
                            species = 'human',
                            ...)
    return(conn)
}
