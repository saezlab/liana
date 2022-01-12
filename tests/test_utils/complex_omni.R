# import the OmniPathR intercell network component
ligrec <- OmnipathR::import_intercell_network(resources =
                                                  c("Cellinker",
                                                    "CellPhoneDB",
                                                    "CellChatDB",
                                                    "ICELLNET",
                                                    "connectomeDB2020",
                                                    "CellTalkDB"
                                                    )
                                              )

ligrec

# Distinct and remove weird ligands
complex_omni <- ligrec %>%
    filter(!category_intercell_source %in% c("activating_cofactor",
                                             "ligand_antagonist")) %>%
    select(source, target, source_genesymbol, target_genesymbol)  %>%
    distinct_at(c("source", "target"), .keep_all = TRUE) %>%
    liana::decomplexify(columns = c("source", "target"))


# False/diplicate interactions alone
xd <- complex_omni %>%
    group_by(source, target) %>%
    mutate(interaction_n = n()) %>%
    filter(interaction_n > 1) %>%
    # which ones are complexes (they are unique)
    mutate(complex_flag = str_detect(source_complex, pattern = "COMPLEX") |
               str_detect(target_complex, pattern = "COMPLEX")) %>%
    filter(!complex_flag)

# We anti join the false interactions
# check those with source, targets > 1
complex_omni2 <- complex_omni %>%
    anti_join(xd) %>%
    group_by(source, target) %>%
    mutate(interaction_n = n())


# TO Filter by loc, topology
xd <- ligrec %>%
    distinct_at(c("source", "target"), .keep_all = TRUE) %>%
    OmnipathR::filter_intercell_network(loc_consensus_percentile = 51,
                                        consensus_percentile = NULL,
                                        transmitter_topology = c('secreted',
                                                                 'plasma_membrane_transmembrane',
                                                                 'plasma_membrane_peripheral'),
                                        receiver_topology = c('plasma_membrane_transmembrane',
                                                              'plasma_membrane_peripheral'),
                                        # min_curation_effort = 1,
                                        ligrecextra = FALSE)
xd$category_intercell_target %>% unique()
