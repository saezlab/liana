#' Function to obtain SingleCellSignalR-like scores
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`LRscore`) as calculated
#'    by SingleCellSignalR
get_sca <- function(lr_res,
                    ...){
    liana_call(
        lr_res = lr_res,
        method = "sca",
        ...
    )
}



#' Function to obtain connectome-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`weight_sc`) as calculated
#'    by Connectome
get_connectome <- function(lr_res,
                           ...){
    liana_call(
        lr_res = lr_res,
        method = "connectome",
        ...
        )
}



#' Function to obtain NATMI-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`edge_specificity`)
#'    as calculated by NATMI
get_natmi <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "natmi",
        ...
    )
}



#' Function to obtain logFC weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with a logFC metric (`logfc_comb`). `logfc_comb` is
#'    calculated as the product of the (1 vs the rest) log2FC for each ligand
#'    and receptor gene
get_logfc <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "logfc",
        ...
    )
}


#' Function to obtain CellPhoneDB-like scores
#' @inheritParams liana_scores
#' @inheritParams cellphonedb_score
#' @inheritDotParams cellphonedb_score
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
#'
#' @details unfinished
get_cellphonedb <- function(lr_res,
                            ...){

    liana_call(
        lr_res = lr_res,
        method = "cellphonedb",
        ...
    )
}



#' Wrapper Function to obtain scores via liana_pipe
#'
#' @inheritParams liana_pipe
#' @inheritDotParams liana_pipe
#' @inheritParams liana_scores
#' @inheritParams recomplexify
#'
#' @export
#'
#' @return lr_res modified to be method-specific
liana_call <- function(method,
                       seurat_object,
                       op_resource,
                       lr_res,
                       decomplexify = TRUE,
                       complex_policy = 'min0',
                       ...){

    liana_scores(.score_specs()[[method]],
                 lr_res = lr_res,
                 complex_policy = complex_policy,
                 decomplexify = decomplexify,
                 ...)
}



#' Function to obtain different scoring schemes
#'
#' @param score_object score_object specific to the test obtained from score_specs
#' @param lr_res ligand-receptor DE results and other stats between clusters
#' @param decomplexify whether to dissociate complexes into subunits and hence
#'    and hence take complexes into account (decomplexify) or not
#' @inheritParams recomplexify
#' @param ... dot params passed to `*_score` functions
#'
#' @return lr_res modified to be method-specific
liana_scores <- function(score_object,
                         lr_res,
                         decomplexify,
                         complex_policy,
                         ...){
    lr_res %<>%
        select(ligand, receptor,
               ends_with("complex"),
               source, target,
               !!score_object@columns)

    if(decomplexify){
        lr_res %<>%
            recomplexify(
                lr_res = .,
                columns = score_object@columns,
                complex_policy = complex_policy)
    }

    args <-
        append(
            list(lr_res = lr_res,
            score_col = score_object@method_score),
            list(...)
        )

    exec(score_object@score_fun, !!!args) %>%
        ungroup()
}




#' Function Used to Calculate the Connectome-like `weight_sc` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#'
#' @noRd
#'
#' @return lr_res with an added `weight_sc` column
connectome_score <- function(lr_res,
                             score_col){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := mean(c(ligand.scaled, receptor.scaled)))
}


#' Function Used to Calculate the NATMI-like `edge_specificity` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @return lr_res with an added `edge_specificity` column
#'
#' @noRd
#'
#' @details In the original NATMI implementation NAs are filtered out, but
#' here replace NAs with 0s, as they are needed to account for complexes
natmi_score <- function(lr_res,
                        score_col){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := ((ligand.expr*(ligand.sum^-1))) *
                   ((receptor.expr*(receptor.sum^-1)))) %>%
        mutate( {{ score_col }} := tidyr::replace_na(.data[[score_col]], 0))
}


#' Function Used to Calculate the logFC products (by default)
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @noRd
#'
#' @return lr_res with an added `logfc_comb` column
logfc_score <- function(lr_res,
                        score_col){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := mean(ligand.log2FC, receptor.log2FC))
}



#' Function Used to Calculate the SigneCellSignalR `LRscore` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#' @param score_col name of the score column
#'
#' @noRd
#'
#' @return lr_res with an added `LRscore` column
sca_score <- function(lr_res,
                      score_col){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} :=
                    (ligand.expr^(1/2) * receptor.expr^(1/2))/
                    (global_mean + ligand.expr^(1/2) * receptor.expr^(1/2))
        )
}

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
#' @param score_col name of the score column
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


### To be finished ----
#' Correlation Coefficient For Interactions
#'
#' @param sce SingleCellExperiment Object
#' @param lr_res a tabble with LR results obtained in the process of liana_pipe
#'
#' @noRd
#'
#' @return
corr_score <- function(lr_res,
                       sce){

    # should filter to cell type A and B first? - this way it's not specific
    # should also remove genes that are not in lr_res to save time
    corr_pairs <- scran::correlatePairs(sce) %>%
        as_tibble()
    corr_pairs <- corr_pairs %>%
        select(gene1=gene2,
               gene2=gene1,
               everything()) %>%
        bind_rows(corr_pairs) %>%
        select(gene1,
               gene2,
               rho,
               corr.FDR=FDR)

    corr_score <- lr_res %>%
        left_join(
            corr_pairs,
            by=c("ligand"="gene1",
                 "receptor"="gene2")
        ) %>%
        distinct()
}
