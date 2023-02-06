#' Wrapper function to run `cell2cell_tensor` with LIANA output.
#'
#' @details This function servers as a one-liner wrapper to the tensor factorisation
#' method described in \href{https://www.nature.com/articles/s41467-022-31369-2}{tensor_cell2cell}.
#' We refer the user to the publication and \href{https://earmingol.github.io/cell2cell/tutorials/ASD/01-Tensor-Factorization-ASD/}{tensor_cell2cell tutorial page}
#' made by the authors. Logically, one should cite cell2cell's paper if their
#' method was used via LIANA.
#'
#' @param sce SingleCellExperiment with `context_df_dict` - i.e. liana results
#' per context (sample) stored at `sce@metadata$liana_res`
#'
#' @param context_df_dict Dictionary (named list) containing a dataframe for
#' each context. The dataframe must contain columns containing sender (source) cells,
#' receiver (target) cells, ligands, receptors, and communication scores, separately.
#' Keys are context names and values are dataframes. NULL by default.
#' If not NULL will be used instead of `sce@metadata$liana_res`.
#'
#' @param sender_col Name of the column containing the sender cells in all context
#' dataframes.
#'
#' @param receiver_col Name of the column containing the receiver cells in all context
#' dataframes.
#'
#' @param ligand_col Name of the column containing the ligands in all context
#' dataframes.
#'
#' @param receptor_col Name of the column containing the receptors in all context
#' dataframes.
#'
#' @param score_col Name of the column containing the communication scores in
#' all context dataframes.
#'
#' @param how  Approach to consider cell types and genes present across multiple contexts.
#' - 'inner' : Considers only cell types and LR pairs that are present in all contexts (intersection).
#' - 'outer' : Considers all cell types and LR pairs that are present across contexts (union).
#' - 'outer_lr' : Considers only cell types that are present in all contexts (intersection),
#'  while all LR pairs that are present across contexts (union).
#' - 'outer_cells' : Considers only LR pairs that are present in all contexts (intersection),
#'  while all cell types that are present across contexts (union).
#'
#' @param lr_fill Value to fill communication scores when a ligand-receptor
#'  pair is not present across all contexts. (NaN by default)
#'
#' @param cell_fill Value to fill communication scores when a cell is not
#'  present across all ligand-receptor pairs or all contexts. (NaN by default)
#'
#' @param lr_sep Separation character to join ligands and receptors into a
#'  LR pair name. ('^' by Default)
#'
#' @param context_order List used to sort the contexts when building the tensor.
#'  Elements must be all elements in `names(context_df_dict)`. (NULL by default)
#'
#' @param sort_elements  Whether alphabetically sorting elements in the
#' InteractionTensor. The Context Dimension is not sorted if a
#' 'context_order' list is provided. (TRUE by default).
#'
#' @param device Device to use when backend is pytorch.
#' Options are: ['cpu', 'cuda:0', None]. NULL by Default
#'
#' @param rank Ranks for the Tensor Factorization (number of factors to deconvolve the original tensor).
#'  If NULL, then rank selection is performed using the `elbow_rank_selection` function.
#'
#' @param seed Random seed integer
#'
#' @param upper_rank Upper bound of ranks to explore with the elbow analysis.
#'
#' @param runs Number of tensor factorization performed for a given rank.
#' Each factorization varies in the seed of initialization. Consider increasing
#' the number of runs, in order to obtain a more robust rank estimate.
#'
#' @param init Initialization method for computing the Tensor Factorization.
#' {‘svd’, ‘random’}
#'
#' @param build_only Whether to only return a tensor instance, without rank
#'   selection and no factorization.
#'
#' @param factors_only whether to return only the factors after factorization
#'
#' @param conda_env name of existing conda environment
#'
#' @param use_available whether to use c2c if available in current env.
#' False by default.
#'
#' @param verbose verbosity logical
#'
#' @param inplace logical (TRUE by default) if liana results are to be saved
#'  to the SingleCellExperiment object (`sce@metadata$liana_res`)
#'
#' @param ... Dictionary containing keyword arguments for the c2c.compute_tensor_factorization function.
#' The function deals with `random_state` (seed) and `rank` internally.
#'
#' @returns an instance of the cell2cell.tensor.BaseTensor class (via reticulate).
#' If build_only is TRUE, then no rank selection or tensor decomposition is returned.
#' Otherwise, returns a tensor with factorization results.
#'
#' @import reticulate basilisk
#'
#' @export
#'
liana_tensor_c2c <- function(sce=NULL,
                             context_df_dict = NULL,
                             score_col = 'LRscore',
                             how='inner',
                             lr_fill=NaN,
                             cell_fill=NaN,
                             lr_sep="^",
                             context_order=NULL,
                             sort_elements=TRUE,
                             device=NULL,
                             rank=NULL,
                             seed = 1337,
                             upper_rank = 25,
                             runs = 3,
                             init = 'svd',
                             build_only = FALSE,
                             factors_only = TRUE,
                             conda_env=NULL,
                             use_available=FALSE,
                             verbose = TRUE,
                             inplace=TRUE,
                             sender_col = "source",
                             receiver_col = "target",
                             ligand_col = "ligand.complex",
                             receptor_col = "receptor.complex",
                             ...){

    if(is.null(context_df_dict)){
        context_df_dict <- sce@metadata$liana_res
    }

    if(is.null(context_df_dict) & is.null(sce)){
        stop("Please provide a valid `sce` OR `context_df_dict` object")
    }

    if(!is.null(context_df_dict) & !is.null(sce)){
        liana_message("`sce` was superseded by `context_df_dict`!",
                      output = "warning",
                      verbose = verbose
                      )
    }

    # Deal with rank
    rank <- if(is.null(rank)){ NULL } else {as.integer(rank)}


    c2c_available <-
    # Load correct conda env
    if(!is.null(conda_env)){
        print(0)
        liana_message(str_glue("Loading `{conda_env}` Conda Environment"),
                      verbose = verbose,
                      output = "message")
        reticulate::use_condaenv(condaenv = conda_env,
                                 conda = "auto",
                                 required = TRUE)

    # check if c2c is available, if not create env with basilisk
    } else if(use_available){
        if(!reticulate::py_module_available("cell2cell")){
            stop("cell2cell was not found")
        }
    } else{
        # load basilisk env
        liana_message(str_glue("Setting up Conda Environment with Basilisk"),
                      verbose = verbose,
                      output = "message")

        # Set up basilisk
        liana_env <- basilisk::BasiliskEnvironment(envname="liana_cell2cell",
                                                   pkgname="liana",
                                                   packages=.lianapy_packages,
                                                   pip=.liana_pips)
        basilisk::basiliskStart(liana_env, testload="scipy.optimize")
    }
    reticulate::py_set_seed(seed)

    # Import c2c
    c2c <- reticulate::import(module = "cell2cell", as="c2c")

    # Format scores to tensor
    liana_message(str_glue("Building the tensor..."),
                  verbose = verbose,
                  output = "message")

    # Build tensor from scores dict
    tensor <- c2c$tensor$dataframes_to_tensor(context_df_dict = context_df_dict,
                                              sender_col = sender_col,
                                              receiver_col = receiver_col,
                                              ligand_col = ligand_col,
                                              receptor_col = receptor_col,
                                              score_col = score_col,
                                              how = how,
                                              lr_fill = lr_fill,
                                              cell_fill = cell_fill,
                                              lr_sep = lr_sep,
                                              context_order = context_order,
                                              order_labels = list('contexts',
                                                                  'interactions',
                                                                  'senders',
                                                                  'receivers'),
                                              sort_elements = sort_elements,
                                              device = device
                                              )

    if(build_only) return(tensor)

    elbow_metric_raw <- NULL

    # estimate factor rank
    if(is.null(rank)){
        liana_message(str_glue("Estimating ranks..."),
                      verbose = verbose,
                      output = "message")
        py$temp <- tensor$elbow_rank_selection(upper_rank=as.integer(upper_rank),
                                               runs=as.integer(runs),
                                               init=init,
                                               automatic_elbow=TRUE,
                                               random_state=as.integer(seed))

        elbow_metric_raw <- tensor$elbow_metric_raw

        rank <- as.integer(tensor$rank)
    }


    # Compute tensor factorization
    liana_message(str_glue("Decomposing the tensor..."),
                  verbose = verbose,
                  output = "message")
    tensor$compute_tensor_factorization(rank = as.integer(rank),
                                        random_state=as.integer(seed),
                                        ...)

    if(factors_only){
        res <- format_c2c_factors(tensor$factors)

        if(!is.null(elbow_metric_raw)){
            res$elbow_metric_raw <- elbow_metric_raw
        }

    } else{
         res <- tensor
        }

    if(!inplace){
        return(res)
    } else{
        sce@metadata$tensor_res <- res
        return(sce)
    }
}


#' Helper function to format factors from c2c output
#'
#' @details The parameters here would correspond to `order_labels`
#'
#' @param factors factors returned from the
#' `cell2cell.tensor.compute_tensor_factorization` function
#'
#' @param contexts name for `contexts` dimension
#'
#' @param interactions name for `interactions` dimension
#'
#' @param senders name for `senders` dimension
#'
#' @param receivers name for `receivers` dimension
#'
#' @param key_sep separator for the sample and condition names. These are
#' used to jointly name the contexts - technically the names of the
#' `context_df_dict` list passed to the `liana_tensor_c2c` function.
#'
#' @noRd
#'
#' @keywords internal
#'
format_c2c_factors <- function(factors,
                               contexts="contexts",
                               interactions="interactions",
                               senders="senders",
                               receivers="receivers"){
    # Format contexts
    factors[[contexts]] %<>%
        set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
        as_tibble(rownames="context") %>%
        mutate(across(where(is.character), ~as.factor(.x)))

    # Format interactions
    factors[[interactions]] %<>%
        set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
        as_tibble(rownames="lr") %>%
        arrange(lr) %>%
        mutate(lr = factor(lr, lr))

    # Format senders
    factors[[senders]] %<>%
        set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
        as_tibble(rownames="celltype") %>%
        mutate(celltype=factor(celltype, celltype))

    # Format receivers
    factors[[receivers]] %<>%
        set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
        as_tibble(rownames="celltype") %>%
        mutate(celltype=factor(celltype, celltype))

    return(factors)
}


#' Returns tensor cell2cell results
#' @param sce SingleCellExperiment with factors output from tensor-cell2cell
#' @param group_col context descriptor - to be obtained from `colData(sce)`
#' @param sample_col context/sample names - obtained from `colData(sce)`
#'
#' @export
get_c2c_factors <- function(sce, group_col=NULL, sample_col){
    factors <- sce@metadata$tensor_res

    if(is.null(group_col)) return(factors)

    factors$contexts <- .join_meta(sce,
                                   factors$contexts,
                                   sample_col = sample_col,
                                   group_col = group_col)
    return(factors)

}


#' Function to plot an Overview of tensor-c2c results
#'
#' @param sce SingleCellExperiment with factors output from tensor-cell2cell
#' @inheritParams get_c2c_factors
#'
#' @export
#'
plot_c2c_overview <- function(sce, group_col, sample_col){

    factors <- get_c2c_factors(sce, group_col, sample_col)

    # Contexts
    contexts <- factors$contexts %>%
        pivot_longer(cols = -c("context", group_col),
                     names_to = "factor", values_to = "loadings"
        ) %>%
        ggplot(aes(x=context, y=loadings, fill=.data[[group_col]])) +
        geom_bar(stat="identity") +
        facet_grid(factor ~ .) +
        theme_bw(base_size = 14) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              strip.text.y = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle('Contexts') +
        ylab(NULL)

    # lr
    lr <- factors$interactions %>%
        pivot_longer(-lr, names_to = "factor", values_to = "loadings") %>%
        ggplot(aes(x=lr, y=loadings)) +
        geom_bar(stat="identity") +
        facet_grid(factor ~ .) +
        theme_bw(base_size = 14) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              strip.background = element_blank(),
              strip.text.y = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle('Interactions') +
        ylab(NULL)


    # Sender cells
    senders <- factors$senders %>%
        pivot_longer(cols = -celltype,
                     names_to = "factor", values_to = "loadings"
        ) %>%
        ggplot(aes(x=celltype, y=loadings,
                   fill=celltype)) +
        geom_bar(stat="identity") +
        facet_grid(factor ~ .) +
        theme_bw(base_size = 14) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              strip.background = element_blank(),
              strip.text.y = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ylab(NULL) +
        ggtitle('Senders')

    # Receiver cells
    receivers <- factors$receivers %>%
        pivot_longer(cols = -celltype,
                     names_to = "factor", values_to = "loadings"
        ) %>%
        ggplot(aes(x=celltype, y=loadings,
                   fill=celltype)) +
        geom_bar(stat="identity") +
        facet_grid(factor ~ .) +
        theme_bw(base_size = 14) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              strip.background = element_blank(),
              axis.ticks.x=element_blank(),
              strip.text.y = element_text(size=15, face = "bold"),
              plot.title = element_text(hjust = 0.5)) +
        ylab(NULL) +
        ggtitle('Receivers')

    # Assemble overview plot
    overview <- patchwork::wrap_plots(list(contexts,
                                           lr,
                                           senders,
                                           receivers
    ),
    ncol=4,
    nrow(1)) +
        patchwork::plot_layout(guides = "collect")

    grid::grid.draw(overview)
}


#' @title Generate boxplots with significance
#'
#' @param sce SingleCellExperiment with factors output from tensor-cell2cell
#' @inheritParams get_c2c_factors
#'
#' @param ... arguments passed to the test used.
#'
#' @export
plot_context_boxplot <- function(sce,
                                 sample_col,
                                 group_col,
                                 test="t.test",
                                 ...){

    factors <- get_c2c_factors(sce,
                               sample_col=sample_col,
                               group_col=group_col)

    ### Alternative
    contexts_data <- factors$contexts %>%
        select(!!all_of(group_col), starts_with("Factor")) %>%
        pivot_longer(-group_col,
                     names_to = "fact",
                     values_to = "loadings") %>%
        mutate(fact = as.factor(fact)) %>%
        group_by(fact) %>%
        group_nest(keep = TRUE) %>%
        mutate(test = map(data,
                            function(.x) broom::tidy(
                                exec(test,
                                     formula=loadings~.x[[group_col]],
                                     data=.x,
                                     ...
                                     ))
                            )
               ) %>%
        unnest(test) %>%
        mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
        mutate(across(where(is.numeric), ~round(.x, digits = 4)))

    all_facts_boxes <- pmap(contexts_data, function(data=data,
                                                    p.adj=p.adj,
                                                    fact=fact,
                                                    ...){
        data %>%
        ggplot(aes(x=.data[[group_col]], y=loadings, color=.data[[group_col]])) +
            geom_boxplot() +
            theme_minimal() +
            ggtitle(fact, str_glue("{str_to_title(test)}, adj.p = {p.adj}"))
    })

    return(all_facts_boxes)

}

#' @title Plot a Heatmap of context loadings
#'
#' @inheritParams plot_context_boxplot
#' @inheritDotParams ComplexHeatmap::Heatmap
#'
#' @returns a ComplexHeatmap object
#'
#' @export
plot_context_heat <- function(sce,
                              sample_col,
                              group_col,
                              ...){

    factors <- get_c2c_factors(sce,
                               sample_col = sample_col,
                               group_col = group_col)

    # Samples dictionary
    meta_dict <- factors$contexts %>%
        select(context, !!group_col) %>%
        deframe()

    contexts_mat <- factors$contexts %>%
        column_to_rownames("context") %>%
        select(starts_with("Factor")) %>%
        t()

    contexts_mat %>%
        ComplexHeatmap::Heatmap(
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                condition = stringr::str_to_title(meta_dict),
                show_annotation_name = FALSE,
                simple_anno_size = grid::unit(0.3, "cm")
                ),
            name = "Context\nloadings",
            ...)

}


#' Function to plot a UMAP of context loadings
#'
#' @inheritParams plot_context_heat
#' @inheritParams dplyr::top_n
#'
#' @inheritDotParams ComplexHeatmap::Heatmap
#'
#' @returns a ComplexHeatmap object
#'
#' @export
#'
plot_lr_heatmap <- function(sce,
                            lr_sep="^",
                            n = 5,
                            ...){

    factors <- get_c2c_factors(sce)

    # Top n Interactions per factor heatmap
    top_lrs <- factors$interactions %>%
        pivot_longer(-lr,
                     names_to = "fact",
                     values_to = "loadings") %>%
        group_by(fact) %>%
        top_n(wt=loadings, n=n) %>%
        pull(lr)

    lrs_mat <- factors$interactions %>%
        filter(lr %in% top_lrs) %>%
        mutate(lr = gsub(as.character(str_glue("\\{lr_sep}")),
                          " -> ", lr)) %>%
        as.data.frame() %>%
        column_to_rownames("lr") %>%
        as.matrix()

    lrs_mat %>%
        ComplexHeatmap::Heatmap(name = "Interaction\nLoadings",
                                ...)
}


#' Generate a geneset resource for each LR
#' @param sce SingleCellExperiment object with liana_tensor_c2c computed
#' @param resource resource with `source`, `target`, `weight` columns
#'
#' @export
#'
#' @returns a tibble in decoupleR format
generate_lr_geneset <- function(sce, resource, lr_sep="^"){

    # Take SingleCellExperiment object as param instead of lrs to avoid the conflict of declaring 'factors' as global variable
    factors <- get_c2c_factors(sce)
    lrs <- factors$interactions %>%
        separate(lr, sep=str_glue("\\{lr_sep}"), into=c("ligand", "receptor"))

    ligand_weights <- assign_lr_weights(lrs,
                                        resource = resource,
                                        entity = "ligand")
    receptor_weights <- assign_lr_weights(lrs,
                                          resource = resource,
                                          entity = "receptor")


    lrs <- lrs %>%
        select(ligand, receptor) %>%
        inner_join(ligand_weights, by="ligand") %>%
        inner_join(receptor_weights, by="receptor") %>%
        rowwise() %>%
        mutate(mor = .set_coh_mean(ligand_weight, receptor_weight,
                                   ligand_source, receptor_source)
        ) %>%
        na.omit() %>%
        dplyr::rename(set=ligand_source) %>%
        select(-receptor_source, -ligand_weight, -receptor_weight) %>%
        ungroup() %>%
        unite(col="lr", sep=lr_sep, ligand, receptor)

    return(lrs)
}




#' check for sign coherency
#' @noRd
#' @keywords internal
.sign_coh_mean <- function(vec){
    if(length(unique(sign(vec))) > 1){
        return(NA)
    } else {
        return(mean(vec))
    }
}

#' Check for sign and gene-set coherency
#' @noRd
#' @keywords internal
.set_coh_mean <- function(x, y, x_set, y_set){
    if(x_set!=y_set){
        return(NA)
    } else{
        weight = .sign_coh_mean(c(x, y))
    }
}


#' Function to calculate gini coefficients for source and target loadings
#'
#' @param loadings loadings for dimension of interest ('senders' or 'receivers')
#' formatted by `format_c2c_factors`
#'
#' @keywords internal
#'
#' @export
calculate_gini <- function(loadings){

    loadings %>%
        set_tidy_names() %>%
        as_tibble(rownames="celltype",
                  .name_repair = ~ vctrs::vec_as_names(...,
                                                       repair = "universal",
                                                       quiet = TRUE)) %>%
        pivot_longer(-celltype, names_to = "factor", values_to = "loadings") %>%
        group_by(factor) %>%
        summarise(gini = .gini(loadings))
}


#' Plot the product of loadings between the source and target loadings
#' within a factor
#'
#' @param factors factors as formatted by `format_c2c_factors`
#'
#' @param factor_of_int factor of interest e.g. Factor.8
#'
#' @inheritDotParams ComplexHeatmap::Heatmap
#'
#' @export
plot_c2c_cells <- function(sce,
                           factor_of_int,
                           ...){
    factors <- get_c2c_factors(sce)

    sender <- factors$senders %>%
        select(celltype, sender_loadings = !!factor_of_int)

    receiver <- factors$receivers %>%
        select(celltype, receiver_loadings = !!factor_of_int)

    comm <- expand_grid(sender = sender$celltype,
                        receiver = receiver$celltype) %>%
        left_join(sender, by=c("sender"="celltype")) %>%
        left_join(receiver, by=c("receiver"="celltype")) %>%
        mutate(loadings = sender_loadings * receiver_loadings)

    commat <- comm %>%
        select(sender, receiver, loadings) %>%
        pivot_wider(id_cols = sender,
                    names_from = receiver,
                    values_from = loadings) %>%
        as.data.frame() %>%
        column_to_rownames("sender") %>%
        as.matrix()

    commat %>%
        ComplexHeatmap::Heatmap(row_title = "Sender",
                                row_names_side = "left",
                                column_title = "Receiver",
                                ...)
}

#' Gini Coefficient Estimate
#' @param x vector
#'
#' @details adapted from DescTools::Gini
#' @noRd
.gini <- function(x){
    if (any(is.na(x)) || any(x < 0)) return(NA_real_)

    n <- length(x)
    x <- sort(x)

    res <- 2 * sum(x * 1:n) / (n*sum(x)) - 1 - (1/n)
    res <- n / (n - 1) * res

    return(max(0, res))
}


# Env vars ----
.lianapy_packages <- c("python==3.8.8",
                       "pandas==1.4.2",
                       "rpy2==3.4.5")

.liana_pips <- c("cell2cell==0.6.3",
                 "anndata2ri==1.0.6")
