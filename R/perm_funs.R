#' Helper Function to generate shuffled means
#'
#' @param lr_res liana_pipe results
#' @param sce SingleCellExperiment Object
#' @param nperms number of permutations
#' @param seed number used to set random seed
#' @inheritParams liana_pipe
#' @inheritParams map_custom
#' @param verbose logical for verbosity
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
                             parallelize = FALSE,
                             workers = 4,
                             assay.type = "logcounts",
                             verbose = TRUE){
    # remove genes absent in lr_res
    lr_genes <- union(lr_res$ligand, lr_res$receptor)
    sce <- sce[rownames(sce) %in% lr_genes, ]

    # shuffle columns
    set.seed(seed)
    shuffled_clusts <- map(1:nperms, function(perm){
        colLabels(sce) %>%
            as_tibble(rownames = "cell") %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })

    # progress_bar
    if(verbose){
        progress_bar <- progress_estimated(nperms)
    } else{
        progress_bar <- NULL
    }


    # generate mean permutations
    perm <- map_custom(.x = shuffled_clusts,
                       .f = mean_permute,
                       sce = sce,
                       assay.type = assay.type,
                       pb = progress_bar,
                       parallelize = parallelize,
                       workers = workers)

    return(perm)
}


#' Function to calculate p-values as in CellPhoneDB
#'
#' @inheritParams get_permutations
#' @param lr_res liana pipe results
#' @param perm_means permutations obtained via `get_permutations`
#' @param score_col name of the score column
#' @param ... placeholder
#'
#' @returns lr_res + pvalue and lr.mean
cellphonedb_score <- function(lr_res,
                              perm_means,
                              parallelize,
                              workers,
                              score_col = "pvalue",
                              verbose = TRUE,
                              ...){
    og_res <- lr_res %>%
        select(ligand, receptor, source, target)

    if(verbose){
        progress_bar <- dplyr::progress_estimated(length(perm_means) * 2)
    } else{
        progress_bar <- NULL
    }

    perm_joined <- perm_means %>%
        map_custom(function(pmean){
            og_res %>%
                distinct() %>%
                liana:::join_means(means = pmean,
                                   source_target = "source",
                                   entity = "ligand",
                                   type = "expr",
                                   pb = progress_bar) %>%
                liana:::join_means(means = pmean,
                                   source_target = "target",
                                   entity = "receptor",
                                   type = "expr",
                                   pb = progress_bar) %>%
                replace(is.na(.), 0) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(lr_mean = mean(c(ligand.expr, receptor.expr)))
        }, parallelize = parallelize, workers = workers) %>%
        bind_rows()

    # calculate quantiles using ecdf null_dists
    quantiles <- lr_res %>%
        group_by(source, target, ligand.complex, receptor.complex) %>%
        mutate(og_mean = mean(c(ligand.expr, receptor.expr))) %>%
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
        mutate(lr.mean = mean(c(ligand.expr, receptor.expr))) %>%
        left_join(pvals_df,
                  by = c("ligand.complex", "receptor.complex",
                         "source", "target")) %>%
        mutate({{ score_col }} := # replace pval of non-expressed rec and ligs
                   ifelse(ligand.expr == 0 || receptor.expr == 0,
                          1,
                          .data[[score_col]]))

    return(lr_res)
}




#' Function to calculate mean LR expression from shuffled cluster label matrices
#'  as done in CellPhoneDB
#'
#' @param sce_matrix single cell expression matrix (transposed)
#' @param col_labels cluster labels
#' @param pb progress bar object
#' @param assay.type assay type (counts, logcounts, etc)
#'
#' @importFrom dplyr progress_estimated
#'
#' @return Returns a list of means per gene calculated with reshuffled
#'    cluster/cell identity labels
mean_permute <- function(col_labels,
                         sce,
                         pb,
                         assay.type){
    if(!is.null(pb)){
        pb$tick()$print()
    }

    scuttle::summarizeAssayByGroup(sce,
                                   ids=col_labels,
                                   statistics = c("mean"),
                                   assay.type = assay.type)@assays@data$mean
}




#' Helper custom map function
#'
#' @inheritParams purrr::map
#' @param parallelize logical whether to parallelize
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

