
#' Helper custom map function
#'
#' @inheritParams purrr::map
#' @param parallelize whether to parallelize
#' @param workers Number of workers to be used in parallelization
#' @param ... params passed to the called function and to map functions
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
#'
#' @param lr_res liana_pipe results
#' @param sce SingleCellExperiment Object
#' @param nperms number of permutations
#' @param seed number used to set random seed
#' @inheritParams liana_pipe
#' @inheritParams map_custom
#'
#' @return Returns a list of shuffled gene means by cluster
#'
#' @details This function could be made generalizable to any set of genes,
#'   depending on the set (currently lr_res genes) that is used to filter - i.e.
#'   it could be replaced with e.g. genes from TF regulons
get_permutations <- function(lr_res,
                             sce,
                             nperms = 1000,
                             seed = 1234,
                             trim = 0.1,
                             parallelize = FALSE,
                             workers = 4,
                             assay.type = "logcounts"){
    # remove genes absent in lr_res
    lr_genes <- union(lr_res$ligand, lr_res$receptor)
    sce <- sce[rownames(sce) %in% lr_genes, ]
    sce_mat <- t(as.matrix(sce@assays@data[[assay.type]]))

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


#' Function to calculate p-values as in CellPhoneDB
#'
#' @inheritParams get_permutations
#' @param lr_res liana pipe results
#' @param score_col name of the score column
#' @param parallelize whether to parallelize
#' @param workers number of workers
#'
#' @returns lr_res + pvalue and lr.mean
cellphonedb_score <- function(lr_res,
                              perm_means,
                              parallelize,
                              workers,
                              score_col = "pvalue"){
    og_res <- lr_res %>%
        select(ligand, receptor, source, target)

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
        }, parallelize = parallelize, workers = workers) %>%
        bind_rows()

    # calculate quantiles using ecdf null_dists
    quantiles <- lr_res %>%
        group_by(source, target, ligand.complex, receptor.complex) %>%
        mutate(og_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
        select(source, target,
               ligand.complex, ligand,
               receptor, receptor.complex,
               og_mean) %>%
        left_join(perm_joined, by = c("source", "target", "ligand", "receptor")) %>%
        group_by(ligand.complex, receptor.complex, source, target) %>%
        group_split() %>%
        map_dbl(function(interaction){
            null_dist <- ecdf(interaction$lr_mean)
            og_mean <- interaction %>% pull("og_mean") %>% unique()
            null_dist(og_mean) %>% pluck(1)
        })

    pvals_df <- lr_res %>%
        group_by(ligand.complex, receptor.complex, source, target) %>%
        group_keys() %>%
        mutate( {{ score_col }} :=  1 - quantiles)

    lr_res %<>%
        rowwise() %>%
        mutate(lr.mean = mean(c(ligand.trunc, receptor.trunc))) %>%
        left_join(pvals_df,
                  by = c("ligand.complex", "receptor.complex",
                         "source", "target")) %>%
        mutate({{ score_col }} := # replace pval of non-expressed rec and ligs
                   ifelse(ligand.trunc == 0 || receptor.trunc == 0,
                          1,
                          .data[[score_col]]))

    return(lr_res)
}




#' Function to calculate mean LR expression from shuffled cluster label matrices
#'  as done in CellPhoneDB
#'
#' @param sce_matrix single cell expression matrix (transposed)
#' @param col_labels cluster labels
#' @param trim truncate ends of mean
#' @param pb progress bar object
#'
#' @importFrom dplyr progress_estimated
#'
#' @return Returns a list of means per gene calculated with reshuffled
#'    cluster/cell identity labels
mean_permute <- function(col_labels,
                         sce_mat,
                         trim,
                         pb){
    pb$tick()$print()

    stats::aggregate(sce_mat,
                     list(col_labels),
                     FUN=mean,
                     trim=trim) %>%
        tibble::as_tibble() %>%
        dplyr::rename(celltype = Group.1) %>%
        tidyr::pivot_longer(-celltype, names_to = "gene") %>%
        tidyr::pivot_wider(names_from=celltype,
                           id_cols=gene,
                           values_from=value) %>%
        tibble::column_to_rownames("gene")
}

