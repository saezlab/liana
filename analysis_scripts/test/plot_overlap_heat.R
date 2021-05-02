CellChatNC <- omni_resources$CellChatDB %>%
    separate(target_genesymbol, into = c("target_genesymbol1",
                                         "target_genesymbol2",
                                         "target_genesymbol3",
                                         "target_genesymbol4")
    ) %>%
    pivot_longer(cols = c("target_genesymbol1",
                          "target_genesymbol2",
                          "target_genesymbol3",
                          "target_genesymbol4"),
                 values_to = "target_genesymbol",
                 names_to = NULL) %>%
    na.omit()


omni_resources <- omni_resources %>%
    purrr::list_modify("Random" = NULL,
                       "Reshuffled" = NULL,
                       "Default" = NULL,
                       "CellChatDB_NC" = CellChatNC)  %>%
    .[order(names(.))]


prepForUpsetRes <- function(named_list){
    map(names(omni_resources), function(l_name){
        omni_resources[[l_name]] %>%
            select(source_genesymbol, target_genesymbol) %>%
            unite("interaction", source_genesymbol, target_genesymbol, sep="_") %>%
            mutate(!!l_name := 1)
    }) %>% reduce(., full_join, by = "interaction") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
        as.data.frame()
}


resources_binary <- prepForUpsetRes(omni_resources)

# Whole universe
jacc_mat <-
    resources_binary %>%
    select(-interaction) %>%
    t() %>%
    get_simil_dist(.,
               sim_dist = "simil",
               method = "Jaccard",
               diag = TRUE) %>%
    as.matrix()
diag(jacc_mat) <- 1

pheatmap(jacc_mat,
         display_numbers = TRUE,
         number_color = "white",
         silent = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         color = viridis::cividis(n=20),
         fontsize = 15,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA,
         na_col="white")

# Replace 1s with the interactions
tmp <- resources_binary %>%
    as_tibble() %>%
    dplyr::mutate(dplyr::across(!starts_with("interaction"),~ifelse(.==1,interaction,.)))




length(intersect(tmp$OmniPath, tmp$CellChatDB))/length(tmp$CellChatDB)
length(intersect(tmp$OmniPath, tmp$CellChatDB))/length(tmp$CellChatDB[tmp$CellChatDB!=0])

length(intersect(tmp$OmniPath, tmp$CellChatDB_NC))/length(tmp$CellChatDB_NC)
length(intersect(tmp$OmniPath, tmp$CellChatDB_NC))/length(tmp$CellChatDB_NC[tmp$CellChatDB_NC!=0])


length(intersect(tmp$OmniPath, tmp$Baccin2019))/length(tmp$Baccin2019)
length(intersect(tmp$OmniPath, tmp$Baccin2019))/length(tmp$Baccin2019[tmp$Baccin2019!=0])


#
tmp2 <- tmp %>% pivot_longer(-interaction,
                             names_to = "resource",
                             values_to = "interact") %>%
    filter(interact != 0) %>%
    select(resource, interaction = interact) %>%
    group_by(resource) %>%
    group_nest() %>%
    mutate(interaction = data %>% map(function(i) i$interaction))



xd <- seq(length(tmp2$interaction)) %>%
    map(function(i)
        map(length(tmp2$interaction), function(j){
            setdiff(tmp2$interaction[i],
                      tmp2$interaction[j])
        })
    )
