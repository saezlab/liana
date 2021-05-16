ligrec <- compile_ligrec_descr()

ligrec_binary <- ligrec %>%
    map(function(res) pluck(res, "interactions") %>%
            distinct_at(.vars = c("target",
                                  "source"))
    ) %>%
    binarize_resources()


binarize_resourcesx <- function(interaction_list){
    map(names(interaction_list), function(l_name){
        interaction_list[[l_name]] %>%
            select(source, target) %>%
            unite("interaction", source, target, sep="_") %>%
            mutate(!!l_name := 1)
    }) %>% reduce(., full_join, by = "interaction") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
        as.data.frame()
}


# Rework
cats <- list(transmitters = "uniprot",
             receivers = "uniprot",
             interactions = c("source",
                              "target")
             )

binarize_resources <- function(entity_list, cols){
    map(names(entity_list), function(l_name){
        entity_list[[l_name]] %>%
            select(!!cols) %>%
            unite("entity", !!cols, sep="_") %>%
            mutate(!!l_name := 1)
    }) %>% reduce(., full_join, by = "entity") %>%
        mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
        mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
        as.data.frame()
}

# i.e. ligrec_binary new
cats %>%
    map2(names(.),
         function(cols, entity){
             ligrec %>% map(function(res)
                 res[[entity]] %>%
                     select(!!cols)
                 ) %>%
                 binarize_resourcesx(cols)
             }) %>% map2(names(.), function(bindata, entity){
                 jaccheat_save(jacc_pairwise(bindata),
                               figure_path(str_glue("{entity}_jaccard_heat.pdf")),
                               "Jaccard Index")


                 overheat_save(interactions_shared(bindata),
                               figure_path(str_glue("{entity}_shared_heat.pdf")),
                               "% Present")
             })

