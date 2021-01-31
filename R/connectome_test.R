library(Seurat)
library(SeuratData)
# library(devtools)
# BiocManager::install("ComplexHeatmap")
# install_github('msraredon/Connectome', ref = 'master')
# install_github('saezlab/OmnipathR')

library(Connectome)
library(OmnipathR)

# Default DB for connectome is Ramilowski

# Install data and load
pbmc3k <- readRDS("input/pbmc3k_processed.rds")

table(Idents(pbmc3k))
Idents(pbmc3k)


# NOTE THIS STEP
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(pbmc3k)]
pbmc3k <- ScaleData(pbmc3k, features = genes)


# HERE NOTE p.values
pbmc3k.con <- CreateConnectome(pbmc3k,
                              species = 'human', # Only relevant for default LR resource
                              min.cells.per.ident = 75,
                              p.values = F,
                              calculate.DOR = F)
# Filter
pbmc3k.con2 <- FilterConnectome(pbmc3k.con,
                                min.pct = 0.1,
                                min.z = 0.25,
                                max.p = 0.05,
                                remove.na = T)

pbmc3k.con2


# Connectome x OmniPath ----
ramilowski <- import_intercell_network(
    interactions_param = list(resources = 'Ramilowski2015'),
    transmitter_param = list(resources = 'Ramilowski2015'),
    receiver_param = list(resources = 'Ramilowski2015')
) %>%
    select("source_genesymbol" ,"target_genesymbol") %>%
    mutate(mode = "xx") %>%
    arrange(.$source_genesymbol) %>%
    as.data.frame()



# connectome ramilowski
fantom5 <- Connectome::ncomms8866_human %>%
    filter(Pair.Evidence %in% c("literature supported")) %>%
    select(2,4, "mode")



head(ramilowski)
head(fantom5)



# Custom Connectome
pbmc3k.cust <- CreateConnectome(pbmc3k,
                               LR.database = 'custom',
                               min.cells.per.ident = 75,
                               p.values = F,
                               calculate.DOR = F,
                               custom.list = ramilowski)

# Filter
pbmc3k.cust2 <- FilterConnectome(pbmc3k.cust,
                                min.pct = 0.1,
                                min.z = 0.25,
                                max.p = 0.05,
                                remove.na = T)




# Connectome_Omni Pipe ----
omni_list <- list('CellPhoneDB',
                  'Ramilowski2015',
                  'Baccin2019',
                  'LRdb',
                  'Kirouac2010',
                  'ICELLNET',
                  'iTALK',
                  'EMBRACE',
                  'HPMR',
                  'Guide2Pharma',
                  'connectomeDB2020',
                  'talklr',
                  'CellTalkDB',
                  'OmniPath')


setClass("OmniCriteria",
         slots=list(interactions_param="list",
                    transmitter_param="list",
                    receiver_param="list"))


# example object
omni_cr_obj <- methods::new("OmniCriteria",
                            interactions_param=list(resources = 'Ramilowski2015'),
                            transmitter_param=list(resources = 'Ramilowski2015'),
                            receiver_param=list(resources = 'Ramilowski2015'))


# example call
import_intercell_network(
    interactions_param = omni_cr_obj@interactions_param,
    transmitter_param = omni_cr_obj@transmitter_param,
    receiver_param = omni_cr_obj@receiver_param)

# empty object
omni_criteria <- methods::new("OmniCriteria",
                              interactions_param=list(),
                              transmitter_param=list(),
                              receiver_param=list())

# Get a list of objects with different resources
xx_res <- omni_list %>%
    map(function(x){
        print(x)

        if(x!="OmniPath"){
            x_obj = methods::new("OmniCriteria",
                                 interactions_param=list(resource=x),
                                 transmitter_param=list(resource=x),
                                 receiver_param=list(resource=x))

            import_intercell_network(
                interactions_param = x_obj@interactions_param,
                transmitter_param = x_obj@transmitter_param,
                receiver_param = x_obj@receiver_param)
        } else{
            import_intercell_network()
        }
    }) %>%
    setNames(omni_list)



xx_conn <- xx_res %>%
    map(function(x){
        conn <- call_connectome(omni_db=x,
                                seurat_obj = pbmc3k,
                                # optional args passed to createConnectom
                                LR.database = 'custom',
                                min.cells.per.ident = 75,
                                p.values = T,
                                calculate.DOR = F) %>%
            FilterConnectome(.,
                             min.pct = 0.1,
                             min.z = 0.25,
                             max.p = 0.05,
                             remove.na = T)
    }) %>%
    setNames(xx_res)




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

    # scale genes to those available in resource
    connectome.genes <- union(lr_db$source_genesymbol, lr_db$target_genesymbol)
    genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj)]
    seurat_obj <- ScaleData(seurat_obj, features = genes)

    # create connectome
    conn <- CreateConnectome(seurat_obj,
                             custom.list = lr_db,
                             ...)
    return(conn)
}


xx_conn[["OmniPath"]]

names(xx_conn)
