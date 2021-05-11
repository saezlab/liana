get_jacc <- function(sig_list, methods, resources){
    get_binary_df(sig_list) %>%
        dplyr::select(ends_with(resources)) %>%
        select(starts_with(methods)) %>%
        filter(rowMeans(.) != 0) %>%
        t() %>%
        get_simil_dist(sim_dist = "simil", "Jaccard")
}

# JI between iTALK and CellChat for Kirouac and ICELLNET
get_jacc(top_lists$top_500,
         c("iTALK", "CellChat"),
         as.character("Kirouac2010", "ICELLNET"))


methods <- c("iTALK", "CellChat", "Squidpy", "Connectome", "NATMI", "SCA")

# pairwise JI between methods
methods_jacc <- methods %>%
    combn(2) %>%
    t() %>%
    as_tibble() %>%
    unite(c("V1", "V2"), col = "method_combo") %>%
    mutate(methods = str_split(method_combo, "_")) %>%
    mutate(jacc_mat = pmap(list(.x = method_combo,
                                .y = methods),
                           .f = function(.x, .y){
                               list(.x = get_jacc(top_lists$top_500,
                                                  .y,
                                                  as.character(get_lr_resources()))
                                    )
                                    }
                           )
           ) %>%
    unnest(jacc_mat) %>%
    rowwise() %>%
    mutate(jacc_mean = mean(jacc_mat)) %>%
    ungroup() %>%
    arrange(desc(jacc_mean))
methods_jacc

# pairwise JI between Resources
resources_jacc <- as.character(get_lr_resources()) %>%
    combn(2) %>%
    t() %>%
    as_tibble() %>%
    unite(c("V1", "V2"), col = "resource_combo") %>%
    mutate(resources = str_split(resource_combo, "_")) %>%
    mutate(jacc_mat = pmap(list(.x = resource_combo,
                                .y = resources),
                           .f = function(.x, .y){
                               list(.x = get_jacc(top_lists$top_500,
                                                  methods,
                                                  .y)
                               )
                           }
    )
    ) %>%
    unnest(jacc_mat) %>%
    rowwise() %>%
    mutate(jacc_mean = mean(jacc_mat)) %>%
    ungroup() %>%
    arrange(desc(jacc_mean))
resources_jacc


# Same Method different resources
methods %>% map(function(met){
    get_jacc(top_lists$top_500,
             met,
             as.character(as.character(get_lr_resources()))) %>%
        mean()
}) %>% setNames(methods)

