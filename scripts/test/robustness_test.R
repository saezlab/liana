# Load prerequisites
library(tidyverse)
library(Seurat)
library(reticulate)

# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat.rds")
# convert labels to factor (SquidPy)
# breast_cancer@meta.data$lt_id <- as.factor(breast_cancer@meta.data$lt_id)

# Get Omni Resrouces
source("scripts/utils/get_omnipath.R")
omni_resources <- get_omni_resources()

# Get Random DB
# source("scripts/utils/shuffle_omnipath.R")
# op_random <- shuffle_omnipath(omni_resources$OmniPath)

# call squidpy
source("scripts/pipes/squidpy_pipe.R")


# DBs to Bench
db_list <- list("CellChatDB" = omni_resources$CellChatDB,
                "CellTalkDB" = omni_resources$CellTalkDB)

# Run with Full Data
subsampling <- c(1, 0.75, 0.5)
subsampled_res <- bench_robust(subsampling,
             call_squidpyR,
             seurat_object = breast_cancer,
             omni_resources= db_list,
             python_path = "/home/dbdimitrov/anaconda3/bin/python",
             ident = "seurat_clusters")  %>%
    purrr::flatten() %>%
    enframe(value="lr_res") %>%
    mutate(name = str_glue("{name}_subsamp_{rep(subsampling, each = length(db_list))}"))

res <- subsampled_res %>%
    filter(str_detect(name, "CellChatDB"))


# Define Truth
ground <- squidpy_default$CellChatDB %>%
    mutate(truth = if_else(pvalue > 0.05 | is.na(pvalue), 0, 1)) %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    select(interaction, truth)

# negative vastly more than positve - need to downsample
summary(as.factor(ground$truth))


# Prepare for ROC
case1 <- squidpy_sub80$CellChatDB  %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    mutate(sub = 0.8)

case2 <- squidpy_sub40$CellChatDB %>%
    unite(ligand, receptor, source, target, col = "interaction") %>%
    mutate(sub = 0.4)

cases <- list("sub_0.8" = case1,
              "sub_0.4" = case2) %>%
    enframe(name = "subsample", value = "lr_results")


df <- cases %>%
    mutate("roc_prep" = lr_results %>% map(function(case){
        left_join(ground, case) %>%
        mutate(response = case_when(truth == 1 ~ 1,
                                    truth == 0 ~ 0)) %>%
            mutate(predictor = abs((pvalue))) %>%
            filter(!is.na(predictor)) %>%
            select(interaction, response, predictor) %>%
            mutate(response = as.factor(response))
    }))

xd <- df %>% mutate(roc = roc_prep %>%
                  map(function(case) calc_curve(case)))


xd <- xd %>%
    select(subsample, roc) %>%
    unnest(roc)



ggplot(xd, aes(x = 1-specificity,
               y = sensitivity,
               colour = subsample)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")





#
# df <- left_join(ground, cases) %>%
#     mutate(response = case_when(truth == 1 ~ 1,
#                                 truth == 0 ~ 0)) %>%
#     mutate(predictor = abs(log2(pvalue))) %>%
#     filter(!is.na(predictor)) %>%
#     select(interaction, response, predictor)
# df$response = factor(df$response, levels = c(1, 0))
#
# df


# ROC
library(yardstick)
res_col_1 <- "sensitivity"
res_col_2 <- "specificity"
curve_fun = yardstick::roc_curve
auc_fun = yardstick::roc_auc

cn = df %>% filter(.data$response == 0)
cp = df %>% filter(.data$response == 1)

r = df %>%
    curve_fun(.data$response, .data$predictor)
auc = df %>%
    auc_fun(.data$response, .data$predictor)

res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
             {{ res_col_2 }} := r %>% pull(res_col_2),
             th = r$.threshold,
             auc = auc$.estimate,
             n = length(which(res$response == 1)),
             cp = nrow(cp),
             cn = nrow(cn)) %>%
    arrange(!!res_col_1, !!res_col_2)








calc_curve(df)

#'
bc_sub80 <- seurat_subsample(seurat_object = breast_cancer)



cases %>%
    mutate(roc = p)






