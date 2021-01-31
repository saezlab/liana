# Load libs
library(Seurat)
library(SeuratData)
library(Connectome)
library(OmnipathR)


# Load seurat object
pbmc3k <- readRDS("input/pbmc3k_processed.rds")

table(Idents(pbmc3k))
Idents(pbmc3k)


# Function to call connectome with customer databases from OmniPath
# returns a filtered connectome df
call_connectome <- function(omni_db,
                            seurat_obj=pbmc3k,
                            ...){
    # Format db to connectome
    lr_db <- omni_db %>%
        select("source_genesymbol" ,"target_genesymbol") %>%
        mutate(mode = "xxx") %>%
        arrange(.$source_genesymbol) %>%
        as.data.frame()

    # scale genes to ligands and receptors available in the resource
    connectome.genes <- union(lr_db$source_genesymbol, lr_db$target_genesymbol)
    genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj)]
    seurat_obj <- ScaleData(seurat_obj, features = genes)

    # create connectome
    conn <- CreateConnectome(seurat_obj,
                             custom.list = lr_db,
                             ...)
    return(conn)
}


# get connectome results
connectome_results <- omni_resources %>%
    map(function(x){
        conn <- call_connectome(omni_db=x,
                                seurat_obj = pbmc3k,
                                # optional args passed to createConnectom
                                LR.database = 'custom',
                                min.cells.per.ident = 75,
                                p.values = F,
                                calculate.DOR = F) %>%
            FilterConnectome(.,
                             min.pct = 0.1,
                             min.z = 0.25,
                             remove.na = T)
    })  %>%
    setNames(names(omni_resources))


# Save results
saveRDS(connectome_results, "output/connectome_res.RDS")

