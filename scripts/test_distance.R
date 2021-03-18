library(tidyverse)
library(broom)

# Binarized Activities
bc_freq <- sig_list %>%
    get_binary_frequencies()


# Read NES
bc_nes <- read.csv("~/Repos/ligrec_decoupleR/input/sc_bc/breast_cancer_seurat323_NES.csv",
                   header = TRUE, row.names = 1) %>%
    `colnames<-`(str_glue("c{seq(0, 12)}")) %>%
    `rownames<-`(str_glue("c{seq(0, 12)}"))

bc_nes
seq(1, 13) %>%
    map(function(x) seq(1, 13) %>% map(function(y)
        str_glue("{bc_nes[x,y]}")))


bc_nes_vec <- map2(.x=str_glue("c{seq(0, 12)}"), .y=seq(1, 13), .f=function(x1, y1){
    map2(.x=str_glue("c{seq(0, 12)}"), .y=seq(1, 13), .f=function(x2, y2){
        str_glue("{x1}_{x2}:{bc_nes[y1, y2]}")
    })
}) %>%
    enframe() %>%
    unnest(value) %>%
    unnest(value) %>%
    separate(value, into = c("clust_pair", "NES"), sep = ":") %>%
    distinct() %>%
    select(clust_pair, NES) %>%
    mutate(NES = as.double(NES))


# join NES and activities
bc_nes_freq <- bc_freq %>%
    mutate(clust_pair = if_else(clust_pair=="c0_12", "c0_c12", clust_pair)) %>%
    left_join(., bc_nes_vec, by = "clust_pair")


# regression of Activity x NES
bc_nes_regression = bc_nes_freq %>%
    mutate(clust_pair = if_else(clust_pair=="c0_12", "c0_c12", clust_pair)) %>%
    group_by(name) %>%
    do(model = lm(freq ~ NES, data = .)) %>%
    mutate(adjr =  model %>% glance() %>% pluck("adj.r.squared")) %>%
    mutate(pval =  model %>% glance() %>% pluck("p.value")) %>%
    select(name, adjr, pval)  %>%
    separate(name, into = c("Method", "Resource"), convert = TRUE) %>%
    mutate_if(is.character, as.factor)


ggplot(bc_nes_regression, aes(x=adjr, y=pval, colour = Method, shape = Resource)) +
    theme_bw(base_size = 26) +
    geom_point(size=5) +
    scale_color_manual(values=brewer.pal(6, "Dark2")) +
    scale_shape_manual(values=1:nlevels(bc_nes_regression$Resource))



# Ranked Frequncies x NES
rank_frequencies


