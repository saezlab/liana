#' Function to Aggregate CCC Method Results
#'
#' @param liana_res LIANA results
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
#' aggragate of the housekeeping interactions of each `external` LIANA++ method.
#'
#'
#' @inheritDotParams RobustRankAggreg::aggregateRanks
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
#' liana_res <- liana_wrap(testdata)
#' # aggregate results from multiple methods
#' liana_res <- liana_aggregate(liana_res)
liana_aggregate <- function(liana_res,
                            resource = NULL,
                            set_cap = "max",
                            cap = NULL,
                            get_ranks = TRUE,
                            get_agrank = TRUE,
                            .score_mode = .score_specs,
                            ...){

    if(!is_tibble(liana_res[[1]]) && is.null(resource)){
        stop("Please provide provide a name for the resource ",
                 "to be plucked and used to aggregate the method results!")
    } else if(!is_tibble(liana_res[[1]])){
        liana_res %<>% map(function(m_results) m_results %>% pluck(resource))
    }

    # fix external methods which return only ligand/receptor, but not .complex
    if(any(names(liana_res) %in% c('cellchat', 'call_natmi',
                                   'call_connectome', 'call_sca',
                                   'call_italk', 'squidpy'))){
        join_cols <- c("source", "target",
                       "ligand", "receptor")
    } else{
        # Default/internal-only liana runs
        join_cols <- c("source", "target",
                       "ligand", "receptor",
                       "ligand.complex", "receptor.complex")
    }


    cap %<>% `%||%`(.select_cap(liana_res, set_cap))

    liana_mlist <- liana_res %>%
        map2(names(.), function(res, method_name){

            if(is.null(.score_mode()[[method_name]])){
                warning(str_glue("Unknown method name or missing specifics for: {method_name}"))
                return()
            } else{
                message(str_glue("Now aggregating {method_name}"))
            }

            method_score <- .score_mode()[[method_name]]@method_score
            desc_order <- .score_mode()[[method_name]]@descending_order

            .method = sym(as.character(str_glue("{method_name}.{method_score}")))
            .rank_col = sym(as.character(str_glue("{method_name}.rank")))

            res %>%
                # split ligand and receptors for cellchat
                {if(method_name=="cellchat")
                    decomplexify(., columns = c("ligand", "receptor")) %>%
                        dplyr::rename(ligand.complex = ligand_complex,
                                      receptor.complex = receptor_complex) else .} %>%
                top_n(n=if_else(desc_order,
                                cap,
                                -cap),
                      wt=!!sym(method_score)) %>%
                mutate( {{ .rank_col }} := .rank_enh(.data[[method_score]], desc_order)) %>%
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
        {`if`(get_ranks, .liana_consensus(., cap))}

    # Get Robust Ranks
    if(get_agrank){
        liana_aggr <- liana_mlist %>%
            .aggregate_rank(., join_cols = join_cols, ...) %>%
            right_join(., liana_aggr,
                       by = join_cols)
    }

    return(liana_aggr)
}

#' Helper function to execute a function on the vector representing the number
#'  of rows in the results for each method
#'
#'  @inheritParams liana::liana_aggregate
#'  @param fun function to execute
#'
.select_cap <- function(liana_res, fun){
    nums <- liana_res %>% map(function(res) nrow(res)) %>% as.numeric
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
.liana_consensus <- function(liana_aggr, cap){
    liana_aggr %>%
        mutate_at(vars(ends_with(".rank")),
              ~ replace(., is.na(.), cap)) %>% # assign .rank_cap to NA
        mutate(mean_rank = pmap_dbl(select(., ends_with(".rank")),
                                    function(...) mean(c(...)))) %>%
        arrange(mean_rank) %>%
        select(source, target,
               ligand, receptor,
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
        min_rank(desc(vec))
    } else{
        min_rank(vec)
    }
}

#' Robust Aggregate ranks using `RobustRankAggreg`
#'
#' @param liana_mlist liana list with method tibbles
#' @param join_cols columns to be concatenated to create entity to be ranked
#'
#' @inheritDotParams RobustRankAggreg::aggregateRanks
.aggregate_rank <- function(liana_mlist, join_cols, ...){
    message("Aggregating Ranks")

    liana_mlist %>%
        map(function(res){
            # bad practice, but almost unavoidable here...
            res %>%
                unite(join_cols,
                      col = "interaction", sep = "⊎") %>%
                pull("interaction")
        }) %>%
        RobustRankAggreg::aggregateRanks(rmat = RobustRankAggreg::rankMatrix(.),
                                         ...) %>%
        as_tibble() %>%
        rename(aggregate_rank = Score,
               interaction = Name) %>%
        separate(col = "interaction", sep = "⊎",
                 into = join_cols)
}

