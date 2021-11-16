require(tidyverse)
require(liana)
require(magrittr)

ligrec <- compile_ligrec(lr_pipeline = FALSE)
ligrec2 <- reform_omni(ligrec)
saveRDS(ligrec2, "inst/omni_resources.rds")


# exclude complexes from OmniPath CCC
ligrec$OmniPath$interactions %<>%
    filter(!(entity_type_intercell_source == "complex" |
                 entity_type_intercell_target == "complex")) %>%
    # filter any mediators
    filter(!(str_detect(category_intercell_source, "cofactor")) &
               !(str_detect(category_intercell_target, "cofactor")) &
               !(str_detect(category_intercell_source, "ligand_regulator")))%>%
    # remove ambiguous/non-membrane associated receptor-receptor interactions
    filter(!(source %in% c("O75462", "Q13261"))) %>%
    # filter any ion_channel/adp-associated interactions
    filter(parent_intercell_target != "ion_channel") %>%
    # interactions need to be reversed
    mutate(target_genesymbol_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                                   "CLEC2A", "CLEC2B", "CLEC2D"),
                                          source_genesymbol,
                                          target_genesymbol),
           sourge_genesymbol_new = ifelse(target_genesymbol %in% c("FGF2", "FGF23", "ALOX5",
                                                                   "CLEC2A", "CLEC2B", "CLEC2D"),
                                          target_genesymbol,
                                          source_genesymbol)) %>%
    dplyr::select(-c(target_genesymbol_new, sourge_genesymbol_new)) %>%
    distinct_at(.vars=c("source", "target"), .keep_all = TRUE)


# Keep only nodes that are part of the interactions
ligrec$OmniPath$receivers %<>%
    dplyr::filter(genesymbol %in% ligrec$OmniPath$interactions$target_genesymbol) %>%
    distinct_at(.vars="genesymbol", .keep_all = TRUE)
ligrec$OmniPath$transmitters %<>%
    dplyr::filter(genesymbol %in% ligrec$OmniPath$interactions$source_genesymbol) %>%
    distinct_at(.vars="genesymbol", .keep_all = TRUE)




# Fix CPDB ----
omni_cpdb <- ligrec$CellPhoneDB$interactions %>%
# check if any ambigous interactions (wrongly annotated ligands/receptors) exist
    rowwise() %>%
    unite(source, target, col = "interaction", remove = FALSE) %>%
    unite(target, source, col = "interaction2", remove = FALSE) %>%
    # identify duplicates
    mutate(dups = if_else(interaction %in% interaction2 |
                              interaction2 %in% interaction,
                          TRUE,
                          FALSE)) %>%
    # ligands which are targets in OmniPath
    mutate(wrong_transitters = (source %in% ligrec$OmniPath$target)) %>%
    # receptors which are ligands in OmniPath
    mutate(wrong_receivers = (target %in% ligrec$OmniPath$source)) %>%
    # filter duplicates which are wrongly annotated
    filter(!(wrong_transitters & dups)) %>%
    filter(!(wrong_receivers & dups))


# check if any non-membrane bound receptors/ligands are regarded
#  as both receptor and ligand and filter sources/ligands which are receptors in OmniPath
duplicated_cpdb <- omni_cpdb %>%
    filter(dups)

# The rest of receptors, we manually check and fix fix
duplicated_cpdb$target %>% unique()
# here we include (secreted) ligand proteins that fall in the target column and hence
# need to be filtered out
mismatched_trasnmitters <- c("P09917")
duplicated_cpdb %>%
    filter(!(target %in% mismatched_trasnmitters)) %>%
    select(interaction, source_genesymbol, target_genesymbol)


# Appropriately Filter CPDB ----
ligrec$CellPhoneDB$interactions %<>%
    # check if any ambigous interactions (wrongly annotated ligands/receptors) exist
    rowwise() %>%
    unite(source, target, col = "interaction", remove = FALSE) %>%
    unite(target, source, col = "interaction2", remove = FALSE) %>%
    # identify duplicates
    mutate(dups = if_else(interaction %in% interaction2 |
                              interaction2 %in% interaction,
                          TRUE,
                          FALSE)) %>%
    # ligands which are targets in OmniPath
    mutate(wrong_transitters = (source %in% ligrec$OmniPath$interactions$target)) %>%
    # receptors which are ligands in OmniPath
    mutate(wrong_receivers = (target %in% ligrec$OmniPath$interactions$source)) %>%
    # filter duplicates which are wrongly annotated
    filter(!(wrong_transitters & dups)) %>%
    filter(!(wrong_receivers & dups)) %>%
    # mismatched transmitters
    filter(!(target %in% c("P09917"))) %>%
    select(-starts_with("interaction"))


# Fix CellChatDB ----
ligrec$CellPhoneDB$interactions %<>%
    # append missing OG CellChatDB interactions
    bind_rows(get_cellchat_missing()) %>%
    mutate(across(where(is_double), ~replace_na(.x, 1))) %>%
    mutate(across(where(is_integer), ~replace_na(.x, 1))) %>%
    mutate(across(where(is_character), ~replace_na(.x, "placeholder")))




## Fix CellCall ----


ligrec$CellCall





# Temporary solution for checks.
saveRDS(ligrec, "inst/omni_resources.rds")




# Improve OmniPath ----
omni <- select_resource("OmniPath")[[1]]

flipped_interactions <- omni %>%
    filter(source %in% target) %>%
    filter((parent_intercell_source %in% c("ligand", "secreted_enzyme")))

# Get sources in flipped:
flipped_interactions$source %>% unique()
# manually checked all of these, and each was a receptor
receptors_in_ligands = c("O75462", "O95971", "P08887",
                         "P10721", "P16871", "P19438",
                         "P20333", "P22455", "P22607",
                         "P24394", "P25445", "P35916",
                         "P35968", "P42702", "Q13261")

# Get target in flipped
flipped_interactions$target %>% unique()
receptors_in_flipped <- c("P26992", "P40189", "P42702", "Q92956", "P09619",
                          "P15509", "P16144", "P19235", "P21926", "P32927",
                          "P36888", "P60033", "Q16827", "Q92729", "Q9HC73",
                          "P20333", "P19438", "P19022", "P54764", "P31785",
                          "P78552", "P00533", "P05556", "P08648", "P35968",
                          "O14786", "P05106", "P17948", "P33151", "P35916",
                          "P41231", "Q12913", "P14784", "P30530")
# ^ all of these are also receptors, since these interactions are not
# categorized as any membraine-located event -> best to filter them

# How to flip interactions:
# flip these interactions in OmniPath? ----
receptors_in_ligands

flipped_interactions %>% names

xd <- flipped_interactions %>%
    select(source, target, source_genesymbol, target_genesymbol)
xd
reflipped_interactions <- xd %>%
    # create temp names to revert
    rename_at(vars(!tidyselect::contains("target")), ~gsub("source", "target2", .x)) %>%
    rename_at(vars(!tidyselect::contains("target2")), ~gsub("target", "source2", .x)) %>%
    # revert symbol and target info
    rename_at(vars(tidyselect::contains("target2")), ~gsub("target2", "target", .x)) %>%
    rename_at(vars(tidyselect::contains("source2")), ~gsub("source2", "source", .x))

xd
reflipped_interactions %>%
    select(names(xd))

# Omni ++ (filtered) ----
omni <- select_resource("OmniPath")[[1]] %>%
    # filter mediators
    filter(!(str_detect(category_intercell_source, "cofactor")) &
               !(str_detect(category_intercell_target, "cofactor")) &
               !(str_detect(category_intercell_source, "ligand_regulator"))) %>%
    # based on the analysis above we remove ambiguous receptor-receptor interactions
    filter(!(source %in% c("O75462", "Q13261")))
