# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds")) %>%
    Seurat::NormalizeData()
op_resource <- select_resource("CellPhoneDB")[[1]]
# Resource Format
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
entity_genes = union(transmitters$gene,
                     receivers$gene)


## CPDB - this part is called in liana_wrap
lr_res <- liana_pipe(seurat_object = seurat_object,
                     op_resource = op_resource,
                     expr_prop = 0.2,
                     trim = 0.1,
                     assay.type = "logcounts")

sce <- seurat_to_sce(seurat_object = seurat_object,
                     entity_genes,
                     assay="RNA")


# custom map to prevent repetition
map_custom <- function(.x, .f, parallelize, workers, ...){
    if(parallelize){
        future::plan(future::multisession, workers = workers)
        furrr::future_map(.x = .x,
                          .f = .f,
                          .options = furrr::furrr_options(seed = TRUE),
                          ...)

    } else{
        purrr::map(.x = .x,
                   .f = .f,
                   ...)
    }
}





#' Helper Function to generate shuffled means
get_pemutations <- function(lr_res,
                            sce,
                            nperms = 1000,
                            seed = 1234,
                            trim = 0.1,
                            parallelize = FALSE,
                            workers = 4){
    # remove genes absent in lr_res
    lr_genes <- union(lr_res$ligand, lr_res$receptor)
    sce <- sce[rownames(sce) %in% lr_genes, ]
    sce_mat <- t(as.matrix(sce@assays@data$counts))

    # shuffle columns
    set.seed(seed)
    shuffled_clusts <- map(1:nperms, function(perm){
        colLabels(sce) %>%
            as_tibble(rownames = "cell") %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })

    # progress_bar
    progress_bar <- progress_estimated(nperms)

    # generate mean permutations
    perm <- map_custom(.x = shuffled_clusts,
                       .f = mean_permute,
                       sce_mat = sce_mat,
                       trim = trim,
                       pb = progress_bar,
                       parallelize = parallelize,
                       workers = workers)

    return(perm)
}


#' cpdb_score
cpdb_score <- function(lr_res,
                       perm_means,
                       parallelize,
                       workers,
                       score_col = "pvalue"){

    og_res <- lr_res %>%
        rowwise() %>%
        mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
        select(ligand, receptor, source, target, og_mean = lr_mean)

    progress_bar <- dplyr::progress_estimated(length(perm_means) * 2)
    perm_joined <- perm_means %>%
        map_custom(function(pmean){
            og_res %>%
                distinct() %>%
                liana:::join_means(means = pmean,
                                   source_target = "source",
                                   entity = "ligand",
                                   type = "trunc",
                                   pb = progress_bar) %>%
                liana:::join_means(means = pmean,
                                   source_target = "target",
                                   entity = "receptor",
                                   type = "trunc",
                                   pb = progress_bar) %>%
                replace(is.na(.), 0) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
        }, parallelize = parallelize, workers = workers)

    pvals_df <- perm_joined %>%
        bind_rows() %>%
        mutate(lr_mean = na_if(lr_mean, 0)) %>%
        group_by(ligand, receptor, source, target) %>%
        mutate(cc = (sum(og_mean >= lr_mean))) %>%
        dplyr::summarise({{ score_col }} :=
                             1 - (sum(og_mean >= lr_mean)/length(perm_means)))

    og_res %>%
        select(ligand, receptor, source, target, lr.mean = og_mean) %>%
        left_join(pvals_df, by = c("ligand", "receptor", "source", "target"))
}


# Run alg
perm_means <- get_pemutations(lr_res,
                              sce,
                              nperms=10,
                              seed=1234,
                              trim=0.1,
                              parallelize = FALSE,
                              workers=4)


liana_cpdb <- cpdb_score(lr_res = lr_res,
                         perm_means = perm_means,
                         parallelize = FALSE,
                         workers = 4,
                         score_col = "pvalue")












##
og_res <- lr_res %>%
    rowwise() %>%
    mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
    select(ligand, receptor, source, target, og_mean = lr_mean)



progress_bar <- dplyr::progress_estimated(length(perm_means) * 2)
perm_joined <- perm_means %>%
    map_custom(function(pmean){
        og_res %>%
            distinct() %>%
            liana:::join_means(means = pmean,
                               source_target = "source",
                               entity = "ligand",
                               type = "trunc",
                               pb = progress_bar) %>%
            liana:::join_means(means = pmean,
                               source_target = "target",
                               entity = "receptor",
                               type = "trunc",
                               pb = progress_bar) %>%
            replace(is.na(.), 0) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
    }, parallelize = TRUE, workers = 4)



pvals_df <- perm_joined %>%
    bind_rows() %>%
    mutate(lr_mean = na_if(lr_mean, 0)) %>%
    group_by(ligand, receptor, source, target) %>%
    mutate(cc = (sum(og_mean >= lr_mean))) %>%
    dplyr::summarise({{ score_col }} := 1 - (sum(og_mean >= lr_mean)/nperms))


liana_res <- og_res %>%
    select(ligand, receptor, source, target, lr_mean = og_mean) %>%
    left_join(pvals_df, by = c("ligand", "receptor", "source", "target"))





# Filter by genes in lr_res
start.time <- Sys.time()

# remove genes absent in lr_res
lr_genes <- union(lr_res$ligand, lr_res$receptor)
sce <- sce[rownames(sce) %in% lr_genes, ]
sce_mat <- t(as.matrix(sce@assays@data$counts))


# shuffle columns
set.seed(seed)
shuffled_clusts <- map(1:nperms, function(perm){
    colLabels(sce) %>%
        as_tibble(rownames = "cell") %>%
        slice_sample(prop=1, replace = FALSE) %>%
        deframe()
})





# generate mean permutations
start.avg <- Sys.time()
if(parallelize){
    future::plan(future::multisession, workers = workers)
    # make into exec and add progress bar
    perm <- furrr::future_map(.x = shuffled_clusts,
                              .f = mean_permute,
                              sce_mat = sce_mat,
                              trim = trim
    )
} else{
    perm <- map(.x = shuffled_clusts,
                .f = mean_permute,
                sce_mat = sce_mat,
                trim = trim
    )
}
end.avg <- Sys.time()


# This becomes cpdb score
# keep only LR_mean
og_res <- lr_res %>%
    rowwise() %>%
    mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
    select(ligand, receptor, source, target, og_mean = lr_mean)

start.join <- Sys.time()
perm_joined <- perm %>%
    map(function(perm_means){ # parallelize + progress
        og_res %>%
            distinct() %>%
            liana:::join_means(means = perm_means,
                               source_target = "source",
                               entity = "ligand",
                               type = "trunc") %>%
            liana:::join_means(means = perm_means,
                               source_target = "target",
                               entity = "receptor",
                               type = "trunc") %>%
            replace(is.na(.), 0) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
    })
end.join <- Sys.time()

pvals_df <- perm_joined %>%
    bind_rows() %>%
    mutate(lr_mean = na_if(lr_mean, 0)) %>%
    group_by(ligand, receptor, source, target) %>%
    mutate(cc = (sum(og_mean >= lr_mean))) %>%
    dplyr::summarise({{ score_col }} := 1 - (sum(og_mean >= lr_mean)/nperms))


liana_res <- og_res %>%
    select(ligand, receptor, source, target, lr_mean = og_mean) %>%
    left_join(pvals_df, by = c("ligand", "receptor", "source", "target")) %>%
    filter(pvalue <=0.05) %>%
    arrange(across(everything())) %>%
    arrange(pvalue)


end.time <- Sys.time()

# Total time
end.time - start.time

# Time to average
end.avg - start.avg

# Time to join
end.join - start.join







## Compare to the other algorithms ----
squidpy_res <- call_squidpy(seurat_object = seurat_object,
                    op_resource = op_resource,
                    threshold = 0.2,
                    slot = "data") %>%
    filter(pvalue <= 0.05) %>%
    arrange(across(everything())) %>%
    arrange(pvalue)

# conda activate cpdb
seurat_object %<>% Seurat::NormalizeData(normalization.method = "LogNormalize")

write.csv(GetAssayData(seurat_object),
            "../cpdb/input/test_em.csv")

write.csv(Idents(seurat_object) %>%
              enframe(name="barcode", value="annotation"),
          file = "../cpdb/input/metadata.csv",
          row.names = FALSE)


## original cpdb
# conda activate cpdb
# cellphonedb method statistical_analysis input/metadata.csv input/test_em.csv --counts-data hgnc_symbol --threshold 0.2
cpdb_means <- read_delim("~/Repos/cpdb/out/means.txt") %>%
    select(-c(gene_a,
              gene_b, secreted, partner_a, partner_b,
              receptor_a, receptor_b, annotation_strategy,
              is_integrin, id_cp_interaction)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(-interacting_pair) %>%
    mutate(value = as.numeric(value))  %>%
    separate(name, into=c("source", "target"), sep = "\\|")

cpdb_pvalues <- read_delim("~/Repos/cpdb/out/pvalues.txt") %>%
    select(-c(gene_a,
              gene_b, secreted, partner_a, partner_b,
              receptor_a, receptor_b, annotation_strategy,
              is_integrin, id_cp_interaction)) %>%
    pivot_longer(-interacting_pair, values_to = "pvalue") %>%
    separate(name, into=c("source", "target"), sep = "\\|") %>%
    arrange(pvalue)

cpdb_res <- cpdb_means %>% left_join(cpdb_pvalues) %>%
    filter(pvalue <= 0.05) %>%
    separate(interacting_pair, into = c("receptor", "ligand"), sep = "_") %>%
    select(ligand, receptor, everything()) %>%
    arrange(across(everything())) %>%
    arrange(pvalue)
# DT

hist(cpdb_res$pvalue)
hist(liana_res$pvalue)
hist(squidpy_res$pvalue)
