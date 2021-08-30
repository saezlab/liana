library(liana)
library(OmnipathR)
library(dplyr)
library(magrittr)
library(purrr)

## Load Prerequisites
liana_path <- system.file(package = "liana")
testdata <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))


## Obtain Surface LRs only
# filter only PM ligs and regs
ligands <-
    OmnipathR::import_omnipath_intercell(
        causality = 'trans',
        topology = c('pmtm', 'pmp'),
        consensus_percentile = 50,
        loc_consensus_percentile = 30
    ) %>%
    pull(genesymbol) %>%
    unique

receptors <-
    OmnipathR::import_omnipath_intercell(
        causality = 'rec',
        topology = c('pmtm', 'pmp'),
        consensus_percentile = 50,
        loc_consensus_percentile = 30
    ) %>%
    pull(genesymbol) %>%
    unique

# get omni and filter
pm_omni <-  select_resource("OmniPath")[[1]] %>%
    filter(source_genesymbol %in% ligands &
               target_genesymbol %in% receptors)

# Run with tools
squidpy_results <- call_squidpy(testdata,
                                op_resource = pm_omni)

connectome_results <- call_connectome(seurat_object = testdata,
                                      op_resource = pm_omni)

# liana Wrap
liana_test <- liana_wrap(seurat_object = testdata,
                         resource='custom',
                         external_resource = pm_omni,
                         cellchat.params=list(.normalize=TRUE))







#### Omni Function + Complex Omni

# We will first filter OmniPath as described in the LIANA paper:
# i) we only retained interactions with literature references,
# ii) we kept interactions only where the receiver protein was plasma membrane transmembrane
# or peripheral according to at least 30% of the localisation annotations
# iii) we only considered interactions between single proteins (interactions between complexes are also available in OmniPath).
op_omni <- generate_omni(resource = 'OmniPath', # this is just necessary in all the calls
                         loc_consensus_percentile = 30,
                         consensus_percentile = NULL,
                         transmitter_topology = c('secreted',
                                                  'plasma_membrane_transmembrane',
                                                  'plasma_membrane_peripheral'),
                         receiver_topology = c('plasma_membrane_transmembrane',
                                               'plasma_membrane_peripheral'),
                         min_curation_effort = 1,
                         ligrecextra = FALSE,
                         remove_complexes = TRUE
                         )

liana_omni <- select_resource("OmniPath")[[1]]
all_equal(op_omni, liana_omni)


op_omni %>% glimpse()










xd$OmniPath %>% {
    if(remove_complexes)
        filter(., !(entity_type_intercell_source == "complex" |
                        entity_type_intercell_target == "complex"))
    else .

}




OmniPath = list(
    interactions = exec(
        intercell_connections,
        !!!op_ic_quality_param,
        !!!op_ia_quality_param
    )
)




OmniPath$interactions %>%
    filter(!(entity_type_intercell_source == "complex" |
                 entity_type_intercell_target == "complex"))


OmnipathR::import_intercell_network() %>%
    OmnipathR::filter_intercell_network(...)


exec(omnipath_intercell,
     !!!op_ic_quality_param,
     !!!op_ia_quality_param)






### Denes suggestions
library(liana)
library(OmnipathR)
library(dplyr)
library(magrittr)
library(purrr)

ligands <-
    OmnipathR::import_omnipath_intercell(
        parent = c('ligand', 'cell_surface_ligand'),
        topology = c('pmtm', 'pmp'),
        consensus_percentile = 50,
        loc_consensus_percentile = 30,
        entity_types = 'protein'
    ) %>%
    pull(genesymbol) %>%
    unique

receptors <-
    OmnipathR::import_omnipath_intercell(
        parent = 'receptor',
        topology = c('pmtm', 'pmp'),
        consensus_percentile = 50,
        loc_consensus_percentile = 30,
        entity_types = 'protein'
    ) %>%
    pull(genesymbol) %>%
    unique

ligand_receptor <-
    liana::compile_ligrec() %>%
    map(
        function(lr){
            lr$transmitters %<>% filter(genesymbol %in% ligands)
            lr$receivers %<>% filter(genesymbol %in% receptors)
            lr$interactions %<>% filter(
                source_genesymbol %in% ligands &
                    target_genesymbol %in% receptors
            )
            return(lr)
        }
    )
