# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)
source("scripts/utils/bench_robust.R")

# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat323.rds")


# Fix for NATMI
clust.anns <- c("c0", "c1","c2",
                "c3","c4","c5",
                "c6", "c7", "c8",
                "c9", "c10", "c11", "c12")
names(clust.anns) <- levels(breast_cancer)
breast_cancer <- RenameIdents(breast_cancer, clust.anns)

# Get Omni Resrouces
# source("scripts/utils/get_omnipath.R")
# omni_resources <- get_omni_resources()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <-readRDS("input/omni_resources.rds")

# Define Subsampling
subsampling <- c(1, 0.8, 0.6, 0.4)


# 1. Squidpy -------------------------------------------------------------------
source("scripts/pipes/squidpy_pipe.R")

# Omni x  Ramilowski x CellPhoneDB
db_list_squidy <- list("OmniPath" = omni_resources$OmniPath,
                       "Ramilowski2015" = omni_resources$Ramilowski2015,
                       "CellPhoneDB" = omni_resources$CellChatDB)

# convert labels to factor (SquidPy)
# breast_cancer@meta.data$lt_id <- as.factor(breast_cancer@meta.data$lt_id)
squidpy_res <- bench_robust(subsampling,
                            call_squidpyR,
                            seurat_object = breast_cancer,
                            omni_resources= db_list_squidy,
                            python_path = "/home/dbdimitrov/anaconda3/bin/python",
                            ident = "seurat_clusters")  %>%
    purrr::flatten() %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("{name}_subsamp_{rep(subsampling, each = length(db_list_squidy))}"))

# saveRDS(squidpy_res, "output/benchmark/squidpy_res.rds")
squidpy_res <- readRDS("output/benchmark/squidpy_res.rds")

squidpy_sub <- squidpy_res %>%
    mutate(name = str_replace(name, "\\_", ",")) %>%
    separate(name, into = c("resource", "subsample"), sep = ",")
squidpy_sub


squidpy_ground <- squidpy_sub %>%
    filter(subsample == "subsamp_1") %>%
    select(resource, lr_res) %>%
    deframe() %>%
    map(function(ground){
        ground %>%
            mutate(truth = if_else(pvalue > 0.05 | is.na(pvalue), 0, 1)) %>%
            unite(ligand, receptor, source, target, col = "interaction") %>%
            select(interaction, truth)
    }) %>%
    enframe(name = "resource", value = "truth")

squidpy_roc <- left_join(squidpy_sub, squidpy_ground, by = "resource") %>%
    rowwise() %>%
    mutate(roc = list(calc_curve(lr_res, truth, predictor_metric = "pvalue")))
# saveRDS(squidpy_roc, "output/benchmark/squidpy_roc.rds")

xd <- squidpy_roc %>%
    unite(resource, subsample, col = "name") %>%
    select(name, roc) %>%
    unnest(roc)


ggplot(xd, aes(x = 1-specificity,
                            y = sensitivity,
                            colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")



# 2. CellChat ------------------------------------------------------------------
source("scripts/pipes/cellchat_pipe.R")


# DBs to Bench
db_list_cellchat <- list("OmniPath" = omni_resources$OmniPath,
                         "Ramilowski2015" = omni_resources$Ramilowski2015,
                         "Default" = NULL)

cellchat_res <- db_list_cellchat %>%
    map(function(db){
        bench_robust(subsampling,
                     lr_call = call_cellchat,
                     op_resource = db,
                     seurat_object = breast_cancer,
                     exclude_anns = c(), # "ECM-Receptor", "Cell-Cell Contact"
                     nboot = 100,
                     thresh = 1,
                     assay = "SCT")
    })
# saveRDS(cellchat_res, "output/benchmark/cellchat_res.rds")

cellchat_sub <- cellchat_res %>%
    purrr::flatten() %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("{rep(names(db_list_cellchat), each = length(subsampling))}_subsamp_{rep(subsampling, times = length(db_list_cellchat))}"))  %>%
    mutate(name = str_replace(name, "\\_", "xx")) %>%
    separate(name, into = c("resource", "subsample"), sep = "xx")



cellchat_ground <- cellchat_sub %>%
    filter(subsample == "subsamp_1") %>%
    select(resource, lr_res) %>%
    deframe() %>%
    map(function(ground){
        ground %>%
            mutate(truth = if_else(pval > 0.05 | is.na(pval), 0, 1)) %>%
            unite(ligand, receptor, source, target, col = "interaction") %>%
            select(interaction, truth)
    }) %>%
    enframe(name = "resource", value = "truth")


cellchat_roc <- left_join(cellchat_sub, cellchat_ground, by = "resource") %>%
    rowwise() %>%
    mutate(roc = list(calc_curve(lr_res, truth, predictor_metric = "pval")))
# saveRDS(cellchat_roc, "output/benchmark/cellchat_roc.rds")

cellchat_rp <- cellchat_roc %>%
    unite(resource, subsample, col = "name") %>%
    select(name, roc) %>%
    unnest(roc) %>%
    ggplot(., aes(x = 1-specificity,
                      y = sensitivity,
                      colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

cellchat_rp



# 3. NATMI ---------------------------------------------------------------------
# call NATMI
source("scripts/pipes/NATMI_pipe.R")
reticulate::repl_python()
py_set_seed(1004)

# save OmniPath Resource to NATMI format
# omni_to_NATMI(omni_resources,
#               omni_path = "input/omnipath_NATMI")

db_list_natmi <- list("OmniPath" = omni_resources$OmniPath,
                      "Ramilowski2015" = omni_resources$Ramilowski2015)

natmi_res <- bench_robust(subsampling,
                          lr_call = call_natmi,
                          omni_resources= db_list_natmi,
                          seurat_object = breast_cancer,
                          omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                          natmi_path = "~/Repos/NATMI",
                          em_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_em.csv",
                          ann_path = "~/Repos/ligrec_decoupleR/input/natmi_subsample/bc_ann.csv",
                          output_path = "~/Repos/ligrec_decoupleR/output/benchmark/natmi",
                          .write_data = TRUE,
                          .subsampling_pipe = TRUE,
                          .default_run = TRUE)
# saveRDS(natmi_res, "output/benchmark/natmi_res.rds")


natmi_sub <- natmi_res %>%
    map(function(x) x %>%
            purrr::list_modify("lrc2a" = NULL)) %>%
    purrr::flatten() %>%
    enframe(name = "resource", value = "lr_res") %>%
    mutate(subsample = str_glue("subsamp_{rep(subsampling, each = 3)}"))



natmi_ground <- natmi_sub %>%
    filter(subsample == "subsamp_1") %>%
    select(resource, lr_res) %>%
    deframe() %>%
    map(function(ground){
        ground %>%
            mutate(top_ntile = ntile(desc(abs(edge_specificity)), 100)) %>%
            mutate(truth = if_else(top_ntile == 100, 1, 0)) %>%
            unite(ligand, receptor, source, target, col = "interaction") %>%
            select(interaction, truth)
    }) %>%
    enframe(name = "resource", value = "truth")

summary(as.factor(natmi_ground$truth[[1]]$truth))


natmi_roc <- left_join(natmi_sub, natmi_ground, by = "resource") %>%
    rowwise() %>%
    mutate(roc = list(calc_curve(lr_res, truth, predictor_metric = "edge_specificity")))
# saveRDS(natmi_roc, "output/benchmark/natmi_roc.rds")


natmi_rp <- natmi_roc %>%
    unite(resource, subsample, col = "name") %>%
    select(name, roc) %>%
    unnest(roc) %>%
    ggplot(., aes(x = 1-specificity,
                  y = sensitivity,
                  colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

natmi_rp




# auroc_tibble %>%
#     select(.data$statistic, .data$auc, .data$filter_crit,
#            .data$set_name, .data$bench_name) %>%
#     unite("name_lvl", .data$set_name, .data$bench_name, .data$filter_crit) %>%
#     pivot_wider(names_from = .data$name_lvl, values_from = .data$auc) %>%
#     column_to_rownames(var = "statistic")  %>%
#     pheatmap(.,
#              cluster_rows = FALSE,
#              cluster_cols = FALSE,
#              treeheight_col = 0,
#              treeheight_row = 0,
#              display_numbers = TRUE,
#              silent = TRUE)



# 4. Connectome ----------------------------------------------------------------
source("scripts/pipes/connectome_pipe.R")

# Freezes with exec/do.call

# Connectome fix - does not work with Seurat4 object
# Requires RNA instead of Spatial as raw counts
db_list_conn <- list("Default" = NULL,
                     "OmniPath" = omni_resources$OmniPath,
                     "Ramilowski2015" = omni_resources$Ramilowski2015)

seurat_subsets <- map(subsampling, function(ss){
    seurat_object = seurat_subsample(breast_cancer, subsampling = ss)
})


connectome_res <- seurat_subsets %>%
    map(function(seurat_sub){
        call_connectome(seurat_object = seurat_sub,
                        op_resource = NULL,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = FALSE,
                        assay = 'SCT')
    }) %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("subsamp_{rep(subsampling, each = length(db_list))}"))



# connectome_results <- omni_resources %>%
#     map(function(op_resource){
#         conn <- call_connectome(op_resource,
#                                 breast_cancer,
#                                 # optional args passed to createConnectome
#                                 min.cells.per.ident = 1,
#                                 p.values = TRUE,
#                                 calculate.DOR = FALSE,
#                                 assay = 'SCT')
#     })  %>%
#     setNames(names(omni_resources))


# Freezes and hangs for me on Seurat::ScaleData when called from exec/do.call
# connectome_res <- bench_robust(subsampling = subsampling,
#                                lr_call = call_connectome,
#                                op_resource = omni_resources$CellChatDB,
#                                seurat_object = breast_cancer,
#                                min.cells.per.ident = 1,
#                                p.values = FALSE,
#                                calculate.DOR = FALSE,
#                                .format = FALSE,
#                                assay = 'SCT')

glimpse(connectome_sub[1,]$lr_res[[1]])

ground <- connectome_sub[1,]$lr_res[[1]] %>%
    mutate(top_ntile = ntile(weight_sc, 1000)) %>%
    mutate(truth = if_else((p_val_adj.lig <= 0.05 & p_val_adj.rec <= 0.05) &
                               top_ntile > 995 &
                               (percent.source > 0.1 & percent.target > 0.1),
                           1,
                           0)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth) %>%
    na.omit()

summary(as.factor(ground$truth))


# attempt to keep the same proportion of hits as other tools
summary(as.factor(ground$truth))
connectome_res$lr_res[[2]]

# Prepare for ROC
roc_res <- connectome_res %>%
    dplyr::filter(!str_detect(name, "_1")) %>% # keep only subsampled
    mutate(roc = lr_res %>% map(function(.df){
        .df %>%
            na.omit() %>%
            calc_curve(., ground, predictor_metric = "weight_sc") # get roc
    }))


xd <- connectome_res$lr_res[[1]]



ggplot(roc_res %>%
           unnest(roc), aes(x = 1-specificity,
                            y = sensitivity,
                            colour = name)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")






# 5. Combine Results ------
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


natmi_heat
squidpy_heat
cellchat_heat

heat_data <- bind_rows(squidpy_heat, cellchat_heat, natmi_heat)

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


