# Get Resources
lr <- compile_ligrec_descr()

# Get mean and median unique interactions, receptors, and ligands per category
lr_overlap <- lr %>% ligrec_overlap
lr_overlap %>%
    map(function(cat) cat %>%
            group_by(resource, unique) %>%
            summarise(unq_or_not = n()) %>%
            ungroup() %>%
            pivot_wider(id_cols = resource,
                        names_from = unique,
                        values_from = unq_or_not) %>%
            mutate_all(~ replace(., is.na(.), 0)) %>%
            bind_rows(summarise_all(filter(.), ~if(is.numeric(.)) sum(.) else "Total")) %>%
            bind_rows(summarise_all(filter(., !(resource %in% c("OmniPath", "Total"))),
                                    ~if(is.numeric(.)) sum(.) else "Total_excl_Omni")) %>%
            mutate(unq_perc = `TRUE`/( `FALSE` + `TRUE`) * 100) %>%
            bind_rows(summarise_all(filter(., !(resource %in% c("OmniPath", "Total", "Total_excl_Omni"))),
                                    ~if(is.numeric(.)) median(.) else "Median_excl_Omni")) %>%
            bind_rows(summarise_all(filter(., !(resource %in% c("Total"))),
                                    ~if(is.numeric(.)) median(.) else "Median"))
    )

