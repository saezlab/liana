# Seurat
library(Seurat)
library(reticulate)
seurat_object <- readRDS("input/pbmc3k_processed.rds")

# Extract data from Seurat ----
# write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object, assay = "RNA", slot = "data"))) - 1),
#           file = "input/test_em.csv",
#           row.names = T)
# write.csv(Idents(object = seurat_object)  %>%
#               enframe(name="barcode", value="annotation"),
#           file = "input/test_metadata.csv",
#           row.names = F)

# Save OmniPath to appropriate CSVs ----
omni_list %>%
    map(function(x){
        write.csv(omni_resources[[x]] %>%
                                  reform_to_NATMI(),
                              file = str_glue("input/omnipath_NATMI/{x}.csv"),
                              row.names = FALSE)
        })


# reform_to_NATMI <-
reform_to_NATMI <- function(op_resource){
    op_resource <- op_resource %>%
        select("Ligand gene symbol" = source_genesymbol,
               "Receptor gene symbol" = target_genesymbol) %>%
        as.data.frame()

    return(op_resource)
}

# Then Copy them to the NATMI DB location


# Call from R
# Arguments:
#   --interDB INTERDB
#                         lrc2p (default) has literature supported ligand-receptor pairs | lrc2a has putative and literature supported ligand-receptor pairs | the user-supplied interaction database can also be used by calling the name of database file without extension
#   --interSpecies INTERSPECIES
#                         human (default) | mouse | expandp | expanda
#   --emFile EMFILE       the path to the normalised expression matrix file with row names (gene identifiers) and column names (cell-type/single-cell identifiers)
#   --annFile ANNFILE     the path to the metafile in which column one has single-cell identifiers and column two has corresponding cluster IDs (see file 'toy.sc.ann.txt' as an example)
#   --species SPECIES     human (default) | mouse | rat | zebrafish | fruitfly | chimpanzee | dog | monkey | cattle | chicken | frog | mosquito | nematode | thalecress | rice | riceblastfungus | bakeryeast | neurosporacrassa | fissionyeast | eremotheciumgossypii | kluyveromyceslactis, 21 species are supported
#   --idType IDTYPE       symbol (default) | entrez(https://www.ncbi.nlm.nih.gov/gene) | ensembl(https://www.ensembl.org/) | uniprot(https://www.uniprot.org/) | hgnc(https://www.genenames.org/) | mgi(http://www.informatics.jax.org/mgihome/nomen/index.shtml) | custom(gene identifier used in the expression matrix)
#   --coreNum CORENUM     the number of CPU cores used, default is one
#   --out OUT             the path to save the analysis results

# Run Omni x NATMI pipe ----
# Set wd to NATMI script and DBs
setwd("~/Repos/NATMI/")

# append default resources to list of omnipath resources
omni_list <- list(
    'CellChatDB',
    'CellPhoneDB',
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
    'OmniPath'
)

natmi_list <- append(omni_list, list("lrc2p", "lrc2a"))

# load sc
expression_matrix = "~/Repos/ligrec_decoupleR/input/test_em.csv"
annotations = "~/Repos/ligrec_decoupleR/input/test_metadata.csv"

# submit native sys requests
natmi_list %>% map(function(resource){
    system(str_glue("python ExtractEdges.py ",
             "--species human ",
             "--emFile {expression_matrix} ",
             "--annFile {annotations} ",
             "--interDB {resource} ",
             "--coreNum 8 ",
             "--out ~/Repos/ligrec_decoupleR/output/NATMI_test/{resource}",
             sep = " "))
    })

# set dir back to project
etwd("~/Repos/ligrec_decoupleR")


natmi_results <- list.files("./output/NATMI_test/",
           all.files = TRUE,
           recursive = TRUE,
           pattern ="Edges_") %>%
    enframe() %>%
    separate(value, into = c("resource", "file"), remove = FALSE) %>%
    mutate(value =  value %>% map(function(csv) read_csv(str_glue("./output/NATMI_test/{csv}")))) %>%
    select(resource, value)

natmi_results
