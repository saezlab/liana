ligrec$OmniPath %>% purrr::modify_at(c("transmitters", "receivers"))



ligrec$OmniPath$receivers %>%
    filter(genesymbol %in% ligrec$OmniPath$interactions$target_genesymbol) %>%
    distinct_at(.vars="genesymbol", .keep_all = TRUE)

ligrec$OmniPath$transmitters %>%
    filter(genesymbol %in% ligrec$OmniPath$interactions$source_genesymbol)

ligrec$OmniPath$interactions$target_genesymbol
ligrec$OmniPath$interactions$source_genesymbol




