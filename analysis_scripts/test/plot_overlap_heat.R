
ligrec <- compile_ligrec_descr()

omni_resources <- ligrec %>%
    map(function(res) pluck(res, "interactions"))

# CellChat No complex
# CellChatNC <- omni_resources$CellChatDB %>%
#     separate(target_genesymbol, into = c("target_genesymbol1",
#                                          "target_genesymbol2",
#                                          "target_genesymbol3",
#                                          "target_genesymbol4",
#                                          "target_genesymbol5")
#     ) %>%
#     pivot_longer(cols = c("target_genesymbol1",
#                           "target_genesymbol2",
#                           "target_genesymbol3",
#                           "target_genesymbol4",
#                           "target_genesymbol5"),
#                  values_to = "target_genesymbol",
#                  names_to = NULL) %>%
#     na.omit()
#
#
# omni_resources <- omni_resources %>%
#     purrr::list_modify("CellChatDB_NC" = CellChatNC) %>%
#     .[order(names(.))]
#



ligrec_binary <- ligrec %>%
    map(function(res) pluck(res, "interactions")) %>%
    binarize_resources()





jacc_mat <-
    ligrec_binary %>%
    select(-interaction) %>%
    t() %>%
    get_simil_dist(.,
               sim_dist = "simil",
               method = "Jaccard",
               diag = TRUE) %>%
    as.matrix()
diag(jacc_mat) <- 1

jacc_df <- jacc_mat %>%
    as.data.frame() %>%
    rownames_to_column("resource")  %>%
    pivot_longer(-resource) %>%
    as.data.frame() %>%
    mutate_at(vars(resource, "name"), list(~recode(., .x=!!!.resource_short)))







interacts_per_resource <- ligrec_binary %>%
    as_tibble() %>%
    dplyr::mutate(dplyr::across(!starts_with("interaction"),
                                ~ifelse(.==1,interaction,.))) %>%
    pivot_longer(-interaction,
                 names_to = "resource",
                 values_to = "interact") %>%
    filter(interact != 0) %>%
    distinct()

intersects_per_resource <- interacts_per_resource %>%
    select(resource, interaction = interact) %>%
    group_by(resource) %>%
    group_nest() %>%
    mutate(interaction = data %>% map(function(i) i$interaction)) %>%
    mutate(intersect = interaction %>% get_intersect(resource)) %>%
    rowwise() %>%
    mutate(resource_len = length(interaction)) %>%
    ungroup() %>%
    unnest(intersect)

shared_per_resource <- intersects_per_resource %>%
    mutate(resource2 = names(intersect)) %>%
    unnest(intersect) %>%
    select(resource, resource2, intersect, resource_len) %>%
    mutate(shared_prop = intersect/resource_len) %>%
    select(resource, resource2, shared_prop) %>%
    pivot_wider(
        id_cols = resource,
        names_from = resource2,
        values_from = shared_prop
    ) %>%
    pivot_longer(-resource) %>%
    as.data.frame()










# Replace 1s with the interactions
tmp <- ligrec_binary %>%
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
    distinct() %>%
    select(resource, interaction = interact) %>%
    group_by(resource) %>%
    group_nest() %>%
    mutate(interaction = data %>% map(function(i) i$interaction)) %>%
    mutate(intersect = interaction %>% get_intersect(resource)) %>%
    rowwise() %>%
    mutate(resource_len = length(interaction)) %>%
    ungroup() %>%
    unnest(intersect) %>%
    mutate(resource2 = names(intersect)) %>%
    unnest(intersect) %>%
    select(resource, resource2, intersect, resource_len)

tmp3 <- tmp2 %>%
    mutate(shared_prop = intersect/resource_len) %>%
    select(resource, resource2, shared_prop) %>%
    pivot_wider(
        id_cols = resource,
        names_from = resource2,
        values_from = shared_prop
        ) %>%
    pivot_longer(-resource) %>%
    as.data.frame()

tmp3

ggplot(data = tmp3, aes(resource, name, fill = value, label = value)) +
        geom_tile() +
    scale_fill_viridis(
        option = 'cividis',
        guide = guide_colorbar(
            title = sprintf('Shared Interactions')
            )
        ) +
    theme_minimal()+
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 16, hjust = 1),
        axis.text.y = element_text(vjust = 1,
                                   size = 16, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
    xlab("Resource") +
    geom_text(aes(resource, name, label = round(value, digits = 3)), color = "white", size = 5)

cairo_pdf("dfsd.pdf", width = 16, height = 9, family = 'DINPro')

print(p)

dev.off()




# --- Diag heatmaps which make no sense


# Get lower triangle of the matrix
get_lower_tri<-function(mat){
    mat[upper.tri(mat)] <- NA
    return(mat)
}

# Get upper triangle of the matrix
get_upper_tri<-function(mat){
    mat[lower.tri(mat)] <- NA
    return(mat)
}


library(reshape2)
tmp3_lowertri <- get_lower_tri(tmp3) %>% melt(na.rm = TRUE)

library(ggplot2)
ggplot(data = tmp3_lowertri, aes(Var1, Var2, fill = value, label = value))+
    geom_tile(color = "white")+
    scale_fill_viridis(
        option = 'cividis',
        guide = guide_colorbar(
            title = sprintf('Shared Interactions')
        )
    ) +
    theme_minimal()+
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 16, hjust = 1),
        axis.text.y = element_text(vjust = 1,
                                   size = 16, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
    xlab("Resource") +
    coord_fixed() +
    geom_text(aes(Var1, Var2, label = round(value, digits = 2)), color = "white", size = 4)




tmp3_uppertri <- get_upper_tri(tmp3) %>% melt(na.rm = TRUE)
ggplot(data = tmp3_uppertri, aes(Var1, Var2, fill = value, label = value))+
    geom_tile(color = "white")+
    scale_fill_viridis(
        option = 'cividis',
        guide = guide_colorbar(
            title = sprintf('Shared Interactions')
        )
    ) +
    theme_minimal()+
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
    xlab("Resource") +
    coord_fixed() +
    geom_text(aes(Var1, Var2, label = round(value, digits = 2)), color = "white", size = 4)

