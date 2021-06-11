require(intercell)
omni_resources <- compile_ligrec(lr_pipeline = TRUE)


#' Helper Function to get uniques values of a given column
#' @param left_res e.g. OmniPath version of a Resource
#' @param right_res e.g. Inbuilt version of a Resource
#' @param column name of the column to check
get_uniques <- function(left_res, right_res, column){
    jacc <- length(intersect(left_res[[column]], right_res[[column]]))/
        length(union(left_res[[column]], right_res[[column]]))
    message(
        str_glue("Jacc Index: {round(jacc, digits = 3)}")
    )

    left_res %>%
        select(!!column) %>%
        filter(!(.data[[column]] %in% right_res[[column]])) %>%
        distinct()
}

# I) Check CellChat
## inbuilt Resource CellChatDB
cc_inbuilt <- CellChat::CellChatDB.human %>%
    pluck("interaction") %>%
    select(interaction_name, ligand, receptor)


cc_omnipath <- cellchat_formatDB(CellChat::CellChatDB.human,
                                 omni_resources$CellChatDB,
                                 exclude_anns = c()) %>%
    pluck("interaction") %>%
    select(interaction_name, ligand, receptor)

# 1,939 interactions for Original DB, 2,551 for OmniPath
cc_omnipath %>% glimpse()
cc_inbuilt %>% glimpse()



# 1) Unique Interactions
# Unique interactions to OmniPath
get_uniques(cc_omnipath, cc_inbuilt, "interaction_name") # Unique to OmniPath
get_uniques(cc_inbuilt, cc_omnipath, "interaction_name") # Unique to Inbuilt

# 2) Unique Ligands
get_uniques(cc_omnipath, cc_inbuilt, "ligand")  # Unique to OmniPath
get_uniques(cc_inbuilt, cc_omnipath, "ligand") # Unique to Inbuilt


# 3) Unique Receptors
get_uniques(cc_omnipath, cc_inbuilt, "receptor") # Unique to OmniPath
get_uniques(cc_inbuilt, cc_omnipath, "receptor") # Unique to Inbuilt



# II) Check SingleCellSignalR
## Inbuilt Resource LRdb
lrdb_omni <- omni_resources$LRdb %>%
    sca_formatDB() %>%
    unite(ligand, receptor, col = "interaction", remove = FALSE)


load("input/LRdb.rda")
lrdb_inbuilt <- LRdb %>%
    unite(ligand, receptor, col = "interaction", remove = FALSE)

lrdb_omni %>% glimpse()
lrdb_inbuilt %>% glimpse()


# 1) Unique Interactions
# Unique interactions to OmniPath
get_uniques(lrdb_omni, lrdb_inbuilt, "interaction") # Unique to OmniPath
get_uniques(lrdb_inbuilt, lrdb_omni, "interaction") # Unique to Inbuilt

# 2) Unique Ligands
get_uniques(lrdb_omni, lrdb_inbuilt, "ligand")  # Unique to OmniPath
get_uniques(lrdb_inbuilt, lrdb_omni, "ligand") # Unique to Inbuilt


# 3) Unique Receptors
get_uniques(lrdb_omni, lrdb_inbuilt, "receptor") # Unique to OmniPath
get_uniques(lrdb_inbuilt, lrdb_omni, "receptor") # Unique to Inbuilt



# III) Check iTALK
## inbuilt Resource iTALK
italk_omni <- omni_resources$iTALK %>%
    italk_formatDB()

italk_inbuilt <- iTALK::database

# Glimpse Both
italk_omni %>% glimpse()
italk_inbuilt %>% glimpse()

# 1) Unique Interactions
# Unique interactions to OmniPath
get_uniques(italk_omni, italk_inbuilt, "Pair.Name") # Unique to OmniPath
get_uniques(italk_inbuilt, italk_omni, "Pair.Name") # Unique to Inbuilt

# 2) Unique Ligands
get_uniques(italk_omni, italk_inbuilt, "Ligand.ApprovedSymbol")  # Unique to OmniPath
get_uniques(italk_inbuilt, italk_omni, "Ligand.ApprovedSymbol") # Unique to Inbuilt


# 3) Unique Receptors
get_uniques(italk_omni, italk_inbuilt, "Receptor.ApprovedSymbol") # Unique to OmniPath
get_uniques(italk_inbuilt, italk_omni, "Receptor.ApprovedSymbol") # Unique to Inbuilt



# IV) Connectome
# Inbuilt Resource is FANTOM5 (aka Ramilowski)
conn_omni <- omni_resources$Ramilowski2015 %>%
    conn_formatDB() %>%
    unite(source_genesymbol, target_genesymbol,
          col = "Pair.Name", remove = FALSE) %>%
    rename("Ligand.ApprovedSymbol" = source_genesymbol,
           "Receptor.ApprovedSymbol" = target_genesymbol)

conn_inbuilt <- Connectome::ncomms8866_human

# Glimpse Both
conn_omni %>% glimpse()
conn_inbuilt %>% glimpse()

# 1) Unique Interactions
# Unique interactions to OmniPath
get_uniques(conn_omni, conn_inbuilt, "Pair.Name") # Unique to OmniPath
get_uniques(conn_inbuilt, conn_omni, "Pair.Name") # Unique to Inbuilt

# 2) Unique Ligands
get_uniques(conn_omni, conn_inbuilt, "Ligand.ApprovedSymbol")  # Unique to OmniPath
get_uniques(conn_inbuilt, conn_omni, "Ligand.ApprovedSymbol") # Unique to Inbuilt


# 3) Unique Receptors
get_uniques(conn_omni, conn_inbuilt, "Receptor.ApprovedSymbol") # Unique to OmniPath
get_uniques(conn_inbuilt, conn_omni, "Receptor.ApprovedSymbol") # Unique to Inbuilt



# V) NATMI
## Inbuilt Resource is ConnectomeDB (aka Ramilowski2015 follow up)
natmi_omni <- read.csv("~/Repos/NATMI/lrdbs/connectomeDB2020.csv") %>%
    unite(Ligand.gene.symbol, Receptor.gene.symbol,
          col = "Pair.Name", remove = FALSE)

natmi_inbuilt <- read.csv("~/Repos/NATMI/lrdbs/lrc2p.csv") %>%
    unite(Ligand.gene.symbol, Receptor.gene.symbol,
          col = "Pair.Name", remove = FALSE)

# Glimpse Both
natmi_omni %>% glimpse()
natmi_inbuilt %>% glimpse()


# 1) Unique Interactions
# Unique interactions to OmniPath
get_uniques(natmi_omni, natmi_inbuilt, "Pair.Name") # Unique to OmniPath
get_uniques(natmi_inbuilt, natmi_omni, "Pair.Name") # Unique to Inbuilt

# 2) Unique Ligands
get_uniques(natmi_omni, natmi_inbuilt, "Ligand.gene.symbol")  # Unique to OmniPath
get_uniques(natmi_inbuilt, natmi_omni, "Ligand.gene.symbol") # Unique to Inbuilt


# 3) Unique Receptors
get_uniques(natmi_omni, natmi_inbuilt, "Receptor.gene.symbol") # Unique to OmniPath
get_uniques(natmi_inbuilt, natmi_omni, "Receptor.gene.symbol") # Unique to Inbuilt
