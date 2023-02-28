#' Function to Aggregate CCC Method Results
#'
#' @param liana_res LIANA results
#' @param aggregate_how way to aggregate, by default (NULL) will aggregate
#'  all passed methods with the approach specified in `liana:::.score_specs`.
#'  Alternative options are `magnitude` and `specificity`.
#' @param set_cap Function used to set ranked cap (i.e. the value that is
#'    assigned to interactions with NA for scores);
#'    By default, this is set to "max", which is the maximum number of interactions
#'    obtained by the methods; Some methods return all possible ligand-receptor
#'    combinations for each possible source and target cell pair - i.e. the
#'    known universe of all possible interactions (based on the CCC resource)
#' @param resource If methods are ran with multiple resources, the name of the
#'    resource of interest needs to be provided
#'    *Note* if a name is not provided, the first results based on the first
#'    resource in the list will be returned
#' @param cap A cap can for all methods can also be manually set, then the top X
#'    interactions, based on the `specificity` scores for each method will be
#'    returned and the ranking will be carried out solely on them
#' @param get_ranks boolean, whether to return consensus ranks for methods
#' @param get_agrank boolean, whether to return aggregate rank using the
#'    `RobustRankAggreg` package.
#' @param .score_mode defines the way that the methods would be aggragate.
#' By default, we use the score of each method which reflects specificity
#' (if available), if not e.g. the case of SCA we use it's sole scoring function.
#' This aggregation is by default done on the basis of the list returns by
#' `.score_mode`. Alternatively, one could pass `.score_housekeep` to obtain an
#' aggragate of the housekeeping interactions of each method.
#' @param join_cols columns by which different method results will be joined.
#' NULL by default, and automatically will handle the columns depending on the
#' methods used.
#' @inheritDotParams .rank_matrix
#'
#' @return Tibble with the interaction results and ranking for each method
#'
#' @details set_cap is the name of the name of a function that is to be executed
#'   on a vector representing the number of rows in the results for each method,
#'   by default this is set to \link{base::max}, but any other function that
#'   works with vectors could be passed - e.g. min, mean, etc.
#'
#' This function also decomplexifies any complex present in the CellChat results
#' which returns complexes by default
#'
#' @export
#'
#' @examples
#' liana_path <- system.file(package = "liana")
#' # load testdata
#' testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
#' # run liana
#' liana_res <- liana_wrap(testdata, method=c("sca", "natmi"))
#' # aggregate results from multiple methods
#' liana_res <- liana_aggregate(liana_res)
liana_aggregate <- function(liana_res,
                            aggregate_how=NULL,
                            resource = NULL,
                            set_cap = "max",
                            cap = NULL,
                            get_ranks = TRUE,
                            get_agrank = TRUE,
                            .score_mode = .score_specs,
                            verbose = TRUE,
                            join_cols = NULL,
                            ...){

    # define approach to aggregate
    if(!is.null(aggregate_how)){
        if(aggregate_how=="magnitude"){
            score_mode = liana:::.score_housekeep()
        } else if(aggregate_how=="specificity"){
            specs <- liana:::.score_specs()
            specs$sca <- NULL # remove SingleCellSignalR score
            specs$call_sca <- NULL # remove SingleCellSignalR score
            score_mode = specs
        } else{
            stop("Please specify an existing aggregate approach!")
        }
    } else{
        score_mode <- .score_mode()
    }

    if(!is_tibble(liana_res[[1]]) && is.null(resource)){
        stop("Please provide provide a name for the resource ",
             "to be plucked and used to aggregate the method results!")
    } else if(!is_tibble(liana_res[[1]])){
        liana_res %<>% map(function(m_results) m_results %>% pluck(resource))
    }

    # fix external methods which return only ligand/receptor, but not .complex
    if(any(startsWith(names(liana_res), "call_"))){
        if(!all(startsWith(names(liana_res), "call_"))){
            liana_message(
                "Using internal and external methods should be done with caution!",
                output = "warning",
                verbose = verbose)
        }

        join_cols %<>% `%||%` (c("source", "target",
                                 "ligand", "receptor"))
    } else{
        # Default/internal-only liana runs
        join_cols %<>% `%||%` (c("source", "target",
                                 "ligand.complex", "receptor.complex"))
    }

    cap %<>% `%||%`(.select_cap(liana_res, set_cap))

    liana_mlist <- liana_res %>%
        map2(names(.), function(res, method_name){

            if(is.null(score_mode[[method_name]])){
                liana_message(
                    str_glue(
                        "Unknown method name or missing specifics for: {method_name}"
                    ), output = "warning",  verbose = verbose)
                return()
            } else{
                liana_message(str_glue("Now aggregating {method_name}"),
                              output = "message",
                              verbose = verbose)
            }

            method_score <- score_mode[[method_name]]@method_score
            desc_order <- score_mode[[method_name]]@descending_order

            .method = sym(as.character(str_glue("{method_name}.{method_score}")))
            .rank_col = sym(as.character(str_glue("{method_name}.rank")))

            res %>%
                # split ligand and receptors for cellchat
                {if(method_name=="call_cellchat")
                    decomplexify(., columns = c("ligand", "receptor")) %>%
                        dplyr::rename(ligand.complex = ligand_complex,
                                      receptor.complex = receptor_complex) else .} %>%
                top_n(n=if_else(desc_order,
                                cap,
                                -cap),
                      wt=!!sym(method_score)) %>%
                mutate( {{ .rank_col }} := .rank_enh(.data[[method_score]],
                                                     desc_order)) %>%
                arrange(!!.rank_col) %>%
                rename( {{ .method }} := method_score) %>%
                select(!!join_cols,
                       !!.method, !!.rank_col) %>%
                mutate(across(c(source, target), as.character)) %>%
                distinct() %>%
                as_tibble()
        }) %>% compact()

    liana_aggr <- liana_mlist %>%
        purrr::reduce(., full_join, by = join_cols) %>%
        {`if`(get_ranks, .liana_consensus(., cap, join_cols))}

    # Get Robust Ranks
    if(get_agrank){
        liana_aggr <- liana_mlist %>%
            .aggregate_rank(join_cols = join_cols,
                            verbose = verbose,
                            ...) %>%
            right_join(., liana_aggr,
                       by = join_cols)
    }

    return(liana_aggr)
}


#' Aggregate CCC Method results and by both magnitude and specificity ranks
#'
#' @param liana_res LIANA results
#' @inheritDotParams liana_aggregate
#'
#' @export
#'
#' @examples
#' liana_path <- system.file(package = "liana")
#' # load testdata
#' testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
#' # run liana
#' liana_res <- liana_wrap(testdata, method=c("sca", "natmi"))
#' # aggregate results from multiple methods
#' liana_res <- rank_aggregate(liana_res)
rank_aggregate <- function(liana_res, ...){
    dots <- list(...)
    if("aggregate_how" %in% names(dots)){
        stop("A consensus for both magnitude and specificity will be assigned!")
    }

    keys <- c("source", "target",
              "ligand.complex",
              "receptor.complex")

    magnitude_rank <- liana_aggregate(liana_res,
                                      aggregate_how="magnitude",
                                      ...) %>%
        dplyr::rename(magnitude_rank = aggregate_rank,
                      magnitude_mean_rank = mean_rank) %>%
        select(!ends_with(".rank"))
    specificity_rank <- liana_aggregate(liana_res,
                                        aggregate_how="specificity",
                                        ...) %>%
        select(all_of(keys),
               specificity_rank = aggregate_rank,
               specificity_mean_rank = mean_rank,
               !ends_with(".rank"))

    consensus_rank <- left_join(magnitude_rank, specificity_rank, by=keys) %>%
        select(all_of(keys), magnitude_rank, specificity_rank, everything()) %>%
        arrange(magnitude_rank)

    return(consensus_rank)
}



#' Helper function to execute a function on the vector representing the number
#'  of rows in the results for each method
#'
#' @param fun function to execute
#' @param liana_res ligand-receptor stats between clusters, output of
#' `liana_pipe`
#'
#' @noRd
.select_cap <- function(liana_res, fun){
    nums <- liana_res %>% map(function(res) nrow(res)) %>% as.numeric()
    exec(fun, nums)
}



#' Get Consensus Rankings for the Methods
#'
#' @param liana_aggr Aggregated method results
#' @param cap Value assigned to NA (by default, the max number of all possible
#'    interactions, depending on the resource)
#'
#' @import purrr tibble
#'
#' @return aggregated liana tibble with consensus ranks
#'
#' @noRd
.liana_consensus <- function(liana_aggr, cap, join_cols){
    liana_aggr %>%
        mutate_at(vars(ends_with(".rank")),
                  ~ replace(., is.na(.), cap)) %>% # assign .rank_cap to NA
        mutate(mean_rank = pmap_dbl(select(., ends_with(".rank")),
                                    function(...) mean(c(...)))) %>%
        arrange(mean_rank) %>%
        select(all_of(join_cols),
               ends_with("_rank"),
               everything())
}


#' Convert to Rank Helper Function
#'
#' @param vec vector to rank
#' @param descending_order boolean, whether to sort in desc
#'
#' @return Rank vector
#' @noRd
.rank_enh <- function(vec, descending_order){
    if(descending_order){
        rank(desc(vec), ties.method = "average")
    } else{
        rank(vec, ties.method = "average")
    }
}

#' Robust Aggregate ranks using `RobustRankAggreg`
#'
#' @param liana_mlist liana list with method tibbles
#' @param join_cols columns to be concatenated to create entity to be ranked
#'
#' @noRd
.aggregate_rank <- function(liana_mlist, join_cols, verbose, ...){
    liana_message("Aggregating Ranks", output = "message", verbose = verbose)

    liana_mlist %>%
        map(function(res){
            # bad practice, but almost unavoidable here...
            res %>%
                unite(join_cols,
                      col = "interaction",
                      sep = "⊎") %>%
                pull("interaction")
        }) %>%
        .rank_matrix %>%
        .robust_rank_agg(.,
                         ...) %>%
        separate(col = "interaction", sep = "⊎",
                 into = join_cols)
}



#' Function to convert a list of characters to a ranked matrix [0,1]
#'
#' @param glist a list of ranked/ordered characters
#'
#' @return a matrix filled with 0-1 values depending on the position order
#' of the characters.
#'
#' @details Adated from Kolde et al., 2012.
#' Required due to the removal of the RobustRankAggregate package from CRAN.
#'
#' @references Kolde, R., Laur, S., Adler, P. and Vilo, J., 2012.
#'  Robust rank aggregation for gene list integration and meta-analysis.
#'  Bioinformatics, 28(4), pp.573-580.
#'
#' @keywords internal
.rank_matrix <- function(glist, ...){

    # Get unique entities
    u.ents <- unique(c(glist, recursive = TRUE))

    # num of entities per col
    num.ents <- length(u.ents)

    # get position for each vect
    pos.mat <- sapply(FUN = match,
                      X = glist,
                      x = u.ents,
                      nomatch = num.ents) %>%
        matrix(nrow = num.ents,
               ncol = length(glist),
               dimnames = list(u.ents, names(glist)))

    # Fill mat /w total ents
    rank.mat <- matrix(num.ents,
                       nrow = num.ents,
                       ncol = length(glist),
                       dimnames = list(u.ents, names(glist)))

    # return rank/total_rank by pos
    return(pos.mat / rank.mat)

}


#' Function to calculate and format aggregate ranks
#'
#' @param rmat ranked matrix formated with `.rank_matrix`
#'
#' @details Adated from Kolde et al., 2012.
#' Required due to the removal of the RobustRankAggregate package from CRAN.
#'
#' @references Kolde, R., Laur, S., Adler, P. and Vilo, J., 2012.
#'  Robust rank aggregation for gene list integration and meta-analysis.
#'  Bioinformatics, 28(4), pp.573-580.
#'
#' @noRd
.robust_rank_agg <- function(rmat){
    tibble(interaction = rownames(rmat), # Name
           # calc aggr ranks
           aggregate_rank = unname(apply(rmat, 1, .rho_scores))) %>% # Score
        arrange(aggregate_rank)
}


#' Calculate (corrected) beta scores and rho values
#'
#' @param r normalized ranks vector [0,1]
#'
#' @return The functions returns a vector of p-values
#'
#' @details Adated from Kolde et al., 2012.
#' Required due to the removal of the RobustRankAggregate package from CRAN.
#'
#' @references Kolde, R., Laur, S., Adler, P. and Vilo, J., 2012.
#'  Robust rank aggregation for gene list integration and meta-analysis.
#'  Bioinformatics, 28(4), pp.573-580.
#'
#' @noRd
.rho_scores <- function(r){
    r <- sort(r)
    n <- length(r) #length is sometimes larger than max -> over-inflates FPs

    # Calc beta p-vals
    p <- pbeta(q=r,
               shape1 = 1:n,
               shape2 = n - (1:n) + 1)

    # correct beta pvals
    .corr_beta_pvals(p = min(p), k=n)
}

#' Correct beta p-vals
#'
#' @param p min p-val
#' @param k number of elements in the vec
#'
#' @noRd
.corr_beta_pvals <- function(p, k){
    min(p * k, 1)
}



#' Helper function to rank each method
#'
#' @param liana_res liana_results for a single method
#' @param method_name name of the method
#' @param mode ranking to be carried out. Accepted modes are `specificity`
#' and `magnitude`. The first is meant to reflect the specificity of interactions
#' across all cell types, while the latter typically reflects how highly expressed
#' is a given interaction.
#'
#' @details this function makes use of liana's `liana:::.score_specs` and
#' `liana:::.score_housekeep` functions.
#'
#' @export
rank_method <- function(liana_res,
                        method_name,
                        mode="specificity"){

    if(mode=="specificity"){
        .score_mode <- liana:::.score_specs
    } else if(mode=="magnitude"){
        .score_mode <- liana:::.score_housekeep
    } else {
        stop("Passed `mode` not found!")
    }

    if(is.null(.score_mode()[[method_name]])){
        stop("Score-method combination not found!")
    }

    method_score <- .score_mode()[[method_name]]@method_score
    desc_order <- .score_mode()[[method_name]]@descending_order

    liana_res %>%
        mutate(rank_col = .rank_enh(.data[[method_score]],
                                    desc_order)) %>%
        arrange(rank_col)
}


