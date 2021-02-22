# 7. Combine Results ------
squidpy_heat <- squidpy_roc %>%
    select(resource, subsample, roc) %>%
    ungroup()  %>%
    mutate(resource = str_replace(resource, "CellPhoneDB", "Default")) %>%
    mutate(subsample = str_replace(subsample, "\\.", ",")) %>%
    mutate(alg = "squidpy") %>%
    unite(col = "key", alg, subsample)
squidpy_heat

natmi_heat <- natmi_roc %>%
    select(resource, subsample, roc) %>%
    ungroup() %>%
    mutate(resource = str_replace(resource, "lrc2p", "Default")) %>%
    mutate(alg = "natmi") %>%
    unite(col = "key", alg, subsample)

cellchat_heat <- cellchat_roc %>%
    select(resource, subsample, roc) %>%
    ungroup() %>%
    mutate(alg = "cellchat") %>%
    unite(col = "key", alg, subsample)

conn_roc

conn_heat <- conn_roc %>%
    mutate(name = str_replace(name, "_", "xx")) %>%
    separate(name, into=c("resource", "subsample"), sep = "xx") %>%
    mutate(subsample = str_replace(subsample, "\\.", ",")) %>%
    mutate(alg = "conn") %>%
    unite(col = "key", alg, subsample) %>%
    select(key, roc, resource)



natmi_heat
squidpy_heat
cellchat_heat
conn_heat

heat_data <- readRDS("output/benchmark/heat_data.rds")
# heat_data <- bind_rows(squidpy_heat, cellchat_heat, natmi_heat)
heat_data <- bind_rows(heat_data, conn_heat)
# saveRDS(heat_data ,"output/benchmark/heat_data.rds")


library(pheatmap)
heatp <- heat_data %>%
    filter(!str_detect(key, "subsamp_1")) %>%
    unnest(roc) %>%
    select(key, auc, resource) %>%
    distinct() %>%
    pivot_wider(names_from = resource, values_from = auc) %>%
    column_to_rownames(var = "key") %>%
    pheatmap(.,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = TRUE,
             silent = TRUE)

heatp


squidpy_roc


