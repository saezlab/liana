#' Wrapper function to run `cell2cell_tensor` with LIANA output.
#'
#' @param context_df_dict A list of tibbles with liana-style results per contexts
#' liana format refers to `liana_aggregate` or `liana_wrap` with 1 method output.
#'
#' @param score_col name of the score to be used (LRscore, by default)
#'
#' @param lr_delim Delimiter that separates ligands and receptors in the
#' ligand-receptor pair in the edgelists (score matrices)
#'
#' @param lr_prop minimum proportions
#'
#' @param lr_fill value to fill in missing interactions across contexts
#' (0 by default)
#' @param lr_min minimum ligand receptor (features) required to keep a
#' cell pair(view)
#'
#' @details This function servers as a one-liner wrapper for factorisation
#' with the method described in \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1}{MOFA}
#'
#' @returns a list with matrices of ligand-receptors (rows) and samples
#' (columns) across cell type pairs (list elements)
#'
#' @export
liana_cc2mofa <- function(context_df_dict,
                          score_col = "LRscore",
                          lr_prop = 0.5,
                          lr_delim = "^",
                          lr_fill = NA,
                          lr_min = 15,
                          cell_delim = "&",
                          lr_cell_delim = ";",
                          key_sep = "|"){


    # Bind all samples scores
    scores <- context_df_dict %>%
        imap(function(res, sample_name) res %>%
                 mutate(sample=sample_name)) %>%
        bind_rows() %>%
        unite(ligand.complex, receptor.complex, sep=lr_delim, col="lr") %>%
        unite(source, target, sep = cell_delim, col="cell_pair") %>%
        select(cell_pair, lr, sample, score_col)


    # Get pairs expressed in at least X prop of samples
    lr <- scores %>%
        pull(lr) %>%
        unique()

    # missing cell_pair-lrs
    cell_pair <- scores %>%
        pull(cell_pair) %>%
        unique()

    # all possible cell_pair -> lr combs
    all_combs <- expand_grid(lr, cell_pair)

    # missing cell_pair-lrs across samples
    missing_pairs <- scores %>%
        split(f = as.factor(.$sample)) %>%
        imap(function(res, sample_name){
            all_combs %>%
                anti_join(res, by = c("cell_pair", "lr")) %>%
                mutate(sample = sample_name)
        }) %>%
        bind_rows()

    # Bind missing pairs, pivot to required format, and replace_na
    scores <- scores %>%
        bind_rows(missing_pairs) %>%
        arrange(sample, lr, cell_pair) %>%
        pivot_wider(names_from = sample,
                    values_from = score_col)
    scores

    scores <- scores %>%
        split(f = as.factor(.$cell_pair)) %>%
        map(function(res){
            view <- res %>%
                mutate(lr_present=rowSums(
                    !is.na(select(., names(context_df_dict)))
                    )) %>%
                filter(lr_present / length(context_df_dict) > lr_prop) %>%
                select(-lr_present) %>%
                mutate(across(names(context_df_dict),
                              ~replace_na(.x, lr_fill))) %>%
                unite(lr, cell_pair, col="key", sep=lr_cell_delim) %>%
                column_to_rownames("key") %>%
                as.matrix()

            if(nrow(view) >= lr_min) view else NULL
        }) %>% compact()

    return(scores)

}

