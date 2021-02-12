# Combine results and plots -------------------------------------------------
# cellchat_results <- readRDS("output/pbmc3k/cellchat_results.rds")
# connectome_results <- readRDS("output/pbmc3k/connectome_results.rds")
# natmi_results <- readRDS("output/pbmc3k/natmi_results.rds")
# squidpy_results <- readRDS("output/pbmc3k/squidpy_results.rds")
omnipath_LRs <- list("cellchat" = cellchat_results$OmniPath,
                     "connectome" = connectome_results$OmniPath,
                     "natmi" = natmi_results$OmniPath,
                     "squidpy" = squidpy_results$OmniPath) %>%
  map(function(x) x %>% filter(!(.$source==.$target))) # filter autocrine



omnipath_LRs$cellchat <- omnipath_LRs$cellchat %>% as_tibble() %>%
  filter(pval < 0.05)


# Already filtered
omnipath_LRs$connectome <- omnipath_LRs$connectome %>% as_tibble()



# I would've used edge_specificity, but this is more similar to what they do
# to compare with CellPhoneDB
# they used an arbitrry threshold to filter
omnipath_LRs$natmi <- omnipath_LRs$natmi  %>%
  mutate(top1p = ntile(edge_avg_expr, 100)) %>%
  filter(top1p > 20) %>%
  filter(edge_specificity > 0.1) %>%
  as_tibble()


omnipath_LRs$squidpy <- omnipath_LRs$squidpy %>%
  filter(pvalue < 0.05) %>%
  as_tibble()



# Default LRs
default_LRs <- list("cellchat_def" = cellchat_results$cellchat_def,
                    "cellchat_omni" = cellchat_results$CellChatDB,
                    "connectome_def" = connectome_results$connectome_def,
                    "connectome_omni" = connectome_results$connectomeDB2020,
                    "natmi_def" = natmi_results$lrc2p, # NATMI uses ConnectomeDB
                    "natmi_omni" = natmi_results$connectomeDB2020,
                    "squidpy" = squidpy_results$CellPhoneDB) %>%
  map(function(x) x %>% filter(!(.$source==.$target)))


# pval
default_LRs$cellchat_def <- default_LRs$cellchat_def %>% as_tibble() %>%
  filter(pval < 0.05)
default_LRs$cellchat_omni <- default_LRs$cellchat_omni %>% as_tibble() %>%
  filter(pval < 0.05)



# already filtered
connectome_results$connectome_def <- connectome_results$connectome_def %>%
  as_tibble()
connectome_results$connectomeDB2020 <- connectome_results$connectomeDB2020 %>%
  as_tibble()


default_LRs$natmi_def <- default_LRs$natmi_def  %>%
  mutate(top1p = ntile(edge_avg_expr, 100)) %>%
  filter(top1p > 20) %>%
  filter(edge_specificity > 0.1) %>%
  as_tibble()

default_LRs$natmi_omni <- default_LRs$natmi_omni %>%
  mutate(top1p = ntile(edge_avg_expr, 100)) %>%
  filter(top1p > 20) %>%
  filter(edge_specificity > 0.1) %>%
  as_tibble()


default_LRs$squidpy <- default_LRs$squidpy %>%
  filter(pvalue < 0.05) %>%
  as_tibble()


# Convert lists to Upset matrix
omni_LRs_matrix <- map(names(omnipath_LRs), function(l_name){
  omnipath_LRs[[l_name]] %>%
    select(1:4) %>%
    unite("interaction", source, target,
          ligand, receptor, sep=".") %>%
    mutate(!!l_name := 1)
}) %>% reduce(., full_join) %>%
  mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
  mutate_at(vars(2:ncol(.)), ~ replace(., . > 0, 1)) %>%
  as.data.frame()

# Upset Plot for each tool with OmniPath
omni_LRs_upset <- upset(omni_LRs_matrix, nsets = ncol(omni_LRs_matrix), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")



default_LRs_matrix <- map(names(default_LRs), function(l_name){
  default_LRs[[l_name]] %>%
    select(1:4) %>%
    unite("interaction", source, target,
          ligand, receptor, sep=".") %>%
    mutate(!!l_name := 1)
}) %>% reduce(., full_join) %>%
  mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
  mutate_at(vars(2:ncol(.)), ~ replace(., . > 0, 1)) %>%
  as.data.frame()

# Upset Plot for each tool with default DB
default_LRs_upset <- upset(default_LRs_matrix, nsets = ncol(default_LRs_matrix), order.by = "freq",
                           point.size = 4, line.size = 2, text.scale	= 2,
                           mainbar.y.label = "Significant Interactions",
                           sets.x.label = "Interactions per tool")
default_LRs_upset

# only default DBs
default_LRs_only <- default_LRs_matrix %>%
  select(interaction, cellchat_def, connectome_def, natmi_def, squidpy)

default_only_upset <- upset(default_LRs_only, nsets = ncol(default_LRs_only), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per tool")
default_only_upset




# 6. Tool Upsets and heatmaps --------------------------------------
## CellChat
cellchat_upset_matrix <- cellchat_results %>%
  map(function(db){
    db %>% filter(pval < 0.05) %>%
      as_tibble()})  %>%
      prepForUpset()
    
    
upset(cellchat_upset_matrix %>% select(-c(cellchat_def)), nsets = ncol(cellchat_upset_matrix), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")

# - Omni
cellchat_upset_plot <- upset(cellchat_upset_matrix %>% select(-c(cellchat_def, OmniPath)), nsets = ncol(cellchat_upset_matrix),
      order.by = "freq", point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")


## NATMI
natmi_upset_matrix <- natmi_results %>%
  map(function(db) db  %>%
        mutate(top1p = ntile(edge_avg_expr, 100)) %>%
        filter(top1p > 20) %>%
        filter(edge_specificity > 0.1) %>%
        as_tibble()) %>%
  prepForUpset()

upset(natmi_upset_matrix %>% select(-c(lrc2a)), nsets = ncol(natmi_upset_matrix), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")

# - Omni
natmi_upset_plot <- upset(natmi_upset_matrix %>% select(-c(lrc2a, lrc2p, OmniPath)), nsets = ncol(natmi_upset_matrix),
      order.by = "freq", point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")



## Connectome
connectome_upset_matrix <- connectome_results %>%
  prepForUpset()

upset(connectome_upset_matrix %>% select(-c(connectome_def)), nsets = ncol(connectome_upset_matrix), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")

# - Omni
connectome_upset_plot <- upset(connectome_upset_matrix %>% select(-c(connectome_def, OmniPath)), nsets = ncol(connectome_upset_matrix),
      order.by = "freq", point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")



## squidpy
squidpy_upset_matrix <- squidpy_results %>%
  map(function(df) df %>% filter(pvalue < 0.05)) %>%
  prepForUpset()

upset(squidpy_upset_matrix, nsets = ncol(squidpy_upset_matrix), order.by = "freq",
      point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")

# - Omni
squidpy_upset_plot <- upset(squidpy_upset_matrix %>% select(-c(OmniPath)), nsets = ncol(squidpy_upset_matrix),
      order.by = "freq", point.size = 4, line.size = 2, text.scale	= 2,
      mainbar.y.label = "Significant Interactions",
      sets.x.label = "Interactions per resource")


# save Upset plot
plot_list <- list("omni_LRs_upset" = omni_LRs_upset,
                  "default_LRs_upset" = default_LRs_upset,
                  "default_only_upset" = default_only_upset,
                  "cellchat_upset_plot" = cellchat_upset_plot,
                  "natmi_upset_plot" = natmi_upset_plot,
                  "connectome_upset_plot" = connectome_upset_plot,
                  "squidpy_upset_plot" = squidpy_upset_plot
                  )



names(plot_list) %>%
  map(function(x){
    png(file= paste(str_glue("output/pbmc3k/{x}.png")),
        width = 1280, height = 920)
    print(plot_list[[x]])
    dev.off()
  })


