#' Wrapper function to run `cell2cell_tensor` with LIANA output.
#'
#' @param lrs_per_sample A list of tibbles with liana-style results per contexts
#' liana format refers to `liana_aggregate` or `liana_wrap` with 1 method output.
#' @param score_col name of the score to be used (LRscore, by default)
#' @param lr_how take the union ("outer") or intersection ("inner") of all
#' interactions across contexts (outer by default)
#' @param cell_how take the union ("outer") or intersection ("inner") of all
#' cell types (outer by default)
#' @param lr_delim Delimiter that separates ligands and receptors in the
#' ligand-receptor pair in the edgelists (score matrices);
#' @param cell_delim Delimiter that separates the communicating source and target
#' cell types in the edgelists (score matrices);
#' @param rank Ranks for the Tensor Factorization (number of factors to
#' deconvolve the original tensor to). If None, then rank selection is performed
#' by estimating the elbow that selects the number of factors according that
#' approximate the full tensor (estimated in terms of normalized error).
#' @param lr_fill value to fill in missing interactions across contexts
#' (0 by default)
#' @param cell_fill value to fill in missing cell types across contexts
#' (0 by default)
#' @param values_fill value to fill missing interactions across cell type pairs
#' within the same context; (0 by default)
#' @param seed Random seed integer
#' @param upper_rank Upper bound of ranks to explore with the elbow analysis
#' @param runs Number of tensor factorization performed for a given rank.
#'  Each factorization varies in the seed of initialization.
#' @param init Initialization method for computing the Tensor Factorization.
#' options include c(‘svd’, ‘random’), 'svd' by default
#' @param build_only whether to simply build the tensor from the contexts
#' and return it as is (i.e. no factorization)
#' @param conda_env name of the conda environment to be used.
#' @param join_cols columns by which we join the contexts
#' @param verbose verbosity boolean
#'
#' @details This function servers as a one-liner wrapper to the tensor factorisation
#' method described in \href{https://www.biorxiv.org/content/10.1101/2021.09.20.461129v2.abstract}{tensor_cell2cell}.
#' We refer the user to the publication and \href{https://earmingol.github.io/cell2cell/tutorials/ASD/01-Tensor-Factorization-ASD/}{tensor_cell2cell tutorial page}
#' made by the authors. Logically, one should cite cell2cell's paper if their
#' method was used via LIANA
#'
#' @returns an instance of the cell2cell.tensor.BaseTensor class (via reticulate).
#' If build_only is TRUE, then no rank selection or tensor decomposition is returned.
#' Otherwise, returns a tensor with factorization results.
#'
#' @export
liana_cc2tensor <- function(lrs_per_sample,
                            score_col = 'LRscore',
                            cell_how = "outer",
                            lr_how = "outer",
                            lr_delim = "^",
                            cell_delim = "⊎",
                            rank = NULL,
                            lr_fill = 0,
                            cell_fill = 0,
                            values_fill = 0,
                            seed = 1337,
                            upper_rank = 50,
                            runs = 1,
                            init = 'svd',
                            build_only = FALSE,
                            scores_only = FALSE,
                            conda_env = "cell2cell",
                            join_cols = c("source", "target",
                                          "ligand.complex", "receptor.complex"),
                            verbose=TRUE){

    # Deal with rank
    rank <- if(is.null(rank)){ NULL } else {as.integer(rank)}

    # Load correct conda env
    liana_message(str_glue("Loading `{conda_env}` Conda Environment"),
                  verbose = verbose,
                  output = "message")
    reticulate::use_condaenv(condaenv = conda_env,
                             conda = "auto",
                             required = TRUE)
    reticulate::py_set_seed(seed)

    reticulate::source_python(
        "~/Repos/liana2/inst/tensor_cell2cell.py" #!!!
        # system.file(package = 'liana', "tensor_cell2cell.py")
    )

    # Convert liana_res across samples/contexts to edgelists
    liana_message(str_glue("Generating and Formatting Score Matrices"),
                  verbose = verbose,
                  output = "message")
    score_matrices <- map(lrs_per_sample,
                          ~.to_edgelist(.x,
                                        join_cols = join_cols,
                                        score_col = score_col,
                                        lr_delim = lr_delim,
                                        cell_delim = cell_delim,
                                        values_fill = values_fill))

    # Columns/Cell Pairs
    cell_pairs <- .format_indices(score_matrices, "colnames", cell_how)
    # Rows/Interactions
    lr_pairs <- .format_indices(score_matrices, "rownames", lr_how)

    # Get celltype order
    cell_order <- str_split(cell_pairs,
                            pattern = cell_delim) %>%
        unlist() %>%
        unique() %>%
        sort()
    cell_order

    # Standardize rows/interactions and cells/columns
    scores <- score_matrices %>%
        # preserve function order!
        map(~.reindex(.x, lr_pairs = lr_pairs,
                      lr_fill = lr_fill) %>%
                .recolumn(., cell_pairs = cell_pairs,
                          cell_fill = cell_fill))

    # Return only scores
    if(scores_only) return(scores)


    # Format scores to tensor and run decomposition
    liana_message(str_glue("Building tensor..."),
                  verbose = verbose,
                  output = "message")
    tensor <- py$get_tensor_cell2cell(scores = scores,
                                      rank = rank,
                                      lr_pairs = lr_pairs,
                                      context_order = names(score_matrices),
                                      cell_order = cell_order,
                                      cell_delim = cell_delim,
                                      lr_delim = lr_delim,
                                      seed = as.integer(seed),
                                      build_only = build_only,
                                      upper_rank = as.integer(upper_rank),
                                      runs = as.integer(runs),
                                      init = init
    )

    return(tensor)
}


#' Helper function to convert a liana_res tibble to an edgelist (score matrix)
#'
#' @param liana_res liana tibble output, typically such from `liana_aggregate`,
#' or `liana_wrap` if using only one method.
#'
#' @param join_cols column by which we join
#' @param score_col the name of the score column that is to be used
#' @param values_fill how to fill any missing ligand receptors (0 by default).
#' The consideration how to fill thse here is that, interaction missing between
#' `cell types` are those which do not occur, hence we fill them with the lowest
#' possible values.
#'
#' @details A score matrix is build by selecting one score from a liana_res
#' object (e.g. aggregate_rank) then widening it via a `pivot_wide`, such
#' that rows represent an interaction, and column cell type pairs.
#'
#' Note that the way that LIANA's recomplexify function works
#' it will return the subunit with the min relevant score and if multiple
#' subunits have the same score, complexes for all of them will be returned.
#' Thus, we join using `ligand.complex` and `receptor.complex`, to remove
#' such redundancies
#'
#' @noRd
.to_edgelist <- function(liana_res,
                         join_cols,
                         score_col,
                         values_fill,
                         lr_delim,
                         cell_delim){
    # Remove redundant subunits
    liana_res %>%
        ### !!! Make a check to see if there are subunits with different scores
        distinct_at(.vars = c(all_of(join_cols), score_col)) %>%
        # to CC matrix (rows are ligand-receptor, columns are cellpairs)
        unite(source, target, col = "cellpairs", sep = cell_delim) %>%
        unite(ligand.complex, receptor.complex, col="lrs", sep = lr_delim) %>%
        pivot_wider(id_cols = lrs,
                    names_from = cellpairs,
                    values_from = !!score_col,
                    values_fill = values_fill) %>%
        column_to_rownames('lrs')
}


#' Helper function to reindex edgelist interactions
#'
#' @details deals with both intersects and unions -
#' determined solely from `lr_pairs`)
#'
#' @param edgelist liana_res formatted by `.to_edgelist`
#' @param lr_pairs ligand-receptor interactions to be considered
#' @param lr_fill policy to deal with missing interactions across the samples
.reindex <- function(edgelist, lr_pairs, lr_fill = NaN){
    edgelist %<>% # filter to only those in intersect
        rownames_to_column('lrs') %>%
        filter(lrs %in% lr_pairs)

    if(!is_empty(setdiff(lr_pairs, edgelist$lrs))){  # i.e. if union

        missing_lrs <- setdiff(lr_pairs, edgelist$lrs) %>%
            enframe(value = 'lrs') %>%
            column_to_rownames('name')

        # Bind missing and fill with 0s
        edgelist <- bind_rows(edgelist, missing_lrs) %>%
            mutate(across(everything(), ~replace_na(.x, lr_fill)))
    }

    return(edgelist %<>%
               arrange(lrs) %>%
               column_to_rownames('lrs'))
}

#' Helper function to format the edgelist columns
#'
#' @details deals with both intersects and unions -
#' determined solely from `cell_pairs`)
#'
#' @param edgelist liana_res formatted by `.to_edgelist`
#' @param cell_pairs cell type pairs to be considered
#' @param cell_fill policy to deal with missing cells across the samples
.recolumn <- function(edgelist, cell_pairs, cell_fill = 0){
    columns = cell_pairs %>%
        map_dfr( ~tibble(!!.x := logical()) )

    diffcols <- setdiff(colnames(columns), colnames(edgelist))

    edgelist %>%
        bind_rows(columns) %>%
        select(all_of(cell_pairs)) %>%
        # Note! I only mutate the missing cols (i.e. cell pairs)
        mutate(across(diffcols, ~replace_na(.x, cell_fill)))
}


#' Helper function to determine how to join the cell types and lrs
#' @noRd
.join_how <- function(how){
    if(how=="inner"){
        intersect
    } else if(how=="outer"){
        union
    } else{
        liana_message(str_glue("{how} is not a valid option."),
                      output = "stop")
    }
}


#' Helper function to format indeces (i.e. lr_pairs and cell_pairs)
#'
#' @param score_matrices `.to_edgelist` output
#' @param what colnames or rownames
#' @param how how to join the indeces across samples (lr_how or cell_how)
.format_indices <- function(score_matrices, what, how){
    what <- get(what)

    idxs <- map(score_matrices, ~what(.x))
    indices <- idxs %>% purrr::reduce(.f = .join_how(how))

    # Preserve order or sort new set (either inner or outer)
    if(all(map_lgl(idxs, ~identical(.x, indices)))){
        indices <- idxs[[1]] # preserve
    } else{
        indices <- unique(indices[order(indices)]) # order
    }

    return(indices)

}

