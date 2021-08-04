#' Function to obtain SingleCellSignalR-like scores
#'
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
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`weight_sc`) as calculated
#'    by connectome
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
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`edge_specificity`)
#'    as calculated by natmi
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
                       decomplexify = TRUE,
                       lr_res,
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
#'
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
            list(...),
            list(lr_res = lr_res,
            score_col = score_object@method_score)
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
#'
#' @noRd
#'
#' @return lr_res with an added `logfc_comb` column
logfc_score <- function(lr_res,
                        score_col){
    lr_res %>%
        rowwise() %>%
        mutate( {{ score_col }} := `*`(ligand.log2FC, receptor.log2FC))
}



#' Function Used to Calculate the SigneCellSignalR `LRscore` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
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


### To be finished ----
#' Function to obtain CellPhoneDB-like scores
#'
#' @inheritDotParams liana_call
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
#'
#' @details unfinished
get_cpdb <- function(seurat_object,
                     op_resource,
                     ...){

    liana_call(
        method = "cpdb",
        seurat_object = seurat_object,
        op_resource = op_resource,
        ...
    )
}



#' Function to calculate pvalues as in CPDB
#'
#' @param lr_res LR_res as returned from `liana_pipe`
#' @param sce single cell object
#'
#' @noRd
#'
#' @return lr_res with pvalues
cpdb_score <- function(lr_res,
                       score_col,
                       sce,
                       nperms = 100,
                       seed = 1234,
                       trim = 0.1,
                       parallelize = FALSE,
                       workers = 4){

    # shuffle columns
    set.seed(seed)
    shuffled_clusts <- map(1:nperms, function(perm){
        colLabels(sce) %>%
            as_tibble(rownames = "cell") %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })

    # generate mean permutations
    if(parallelize){
        future::plan(future::multisession, workers = workers)
        perm <- furrr::future_map(.x = shuffled_clusts,
                                  .f = cpdb_permute,
                                  sce_mat = sce_mat,
                                  trim = trim,
                                  .progress=TRUE
                                  )
    } else{
        perm <- map(.x = shuffled_clusts,
                    .f = cpdb_permute,
                    sce_mat = sce_mat,
                    trim = trim
                    )
    }

    # keep only LR_mean
    lr_og <- lr_res %>%
        rowwise() %>%
        mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
        select(ligand, receptor, source, target, og_mean = lr_mean)

    perm %<>%
        map(function(perm_means){
            lr_og %>%
                liana:::join_means(means = perm_means,
                                   source_target = "source",
                                   entity = "ligand",
                                   type = "trunc") %>%
                liana:::join_means(means = perm_means,
                                   source_target = "target",
                                   entity = "receptor",
                                   type = "trunc") %>%
                dplyr::rowwise() %>%
                dplyr::mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
        })


    pvals_df <- perm %>%
        bind_rows() %>%
        group_by(ligand, receptor, source, target) %>%
        mutate(lr_mean = na_if(lr_mean, 0)) %>%
        summarise({{ score_col }} := 1 - (sum(og_mean >= lr_mean)/nperms))


    lr_res %>%
        select(ligand, receptor, source, target, lr_mean) %>%
        left_join(pvals_df, by = c("ligand", "receptor", "source", "target"))
}




#' Function to calculate mean LR expression from shuffled cluster label matrices
#'  as done in CellPhoneDB
#'
#' @param sce_matrix single cell expression matrix (transposed)
#' @param col_labels cluster labels
#' @param trim truncate ends of mean
#' @param lr_og LR res with only relevant
#'
#' @return Returns a list of means per gene calculated with reshuffled
#'    cluster/cell identity labels
cpdb_permute <- function(col_labels,
                         sce_mat,
                         trim){

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

