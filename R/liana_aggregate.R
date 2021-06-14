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
#'
#' @return Tibble with the interaction results and ranking for each method
#'
#' @details set_cap is the name of the name of a function that is to be executed
#'   on a vector representing the number of rows in the results for each method,
#'   by default this is set to \link{base::max}, but any other function that
#'   works with vectors could be passed - e.g. min, mean, etc.
liana_aggregate <- function(liana_res,
                            resource = NULL,
                            set_cap = "max",
                            cap = NULL){

    if(!is_tibble(liana_res[[1]]) && is.null(resource)){
        stop("Please provide provide a name for the resource, ",
                 "otherwise the first resource will be plucked!")
    } else if(!is_tibble(liana_res[[1]])){
        liana_res %<>% map(function(m_results) m_results %>% pluck(resource))
    }

    cap %<>% `%||%`(.select_cap(liana_res, set_cap))

    liana_res %>%
        map2(names(.), function(res, method_name){
            method_score <- .rank_specs[[method_name]]@method_score
            parm_order <- .rank_specs[[method_name]]@descending_order
            .method = sym(as.character(str_glue("{method_name}.{method_score}")))
            .rank_col = sym(as.character(str_glue("{method_name}.rank")))

            res %>%
                top_n(n=if_else(parm_order,
                                cap,
                                -cap),
                      wt=!!sym(method_score)) %>%
                mutate( {{ .rank_col }} := dense_rank(desc(!!sym(method_score)))) %>%
                arrange(!!.rank_col) %>%
                rename( {{ .method }} := method_score) %>%
                select(source, ligand, target, receptor, !!.method, !!.rank_col) %>%
                distinct() %>%
                as_tibble()
        }) %>%
        purrr::reduce(., full_join, by = c("source", "ligand", # Full join all results
                                           "target", "receptor"))
}

#' Helper function to execute a function on the vector representing the number
#'  of rows in the results for each method
#'  @inheritParams liana::liana_aggregate
#'  @param fun function to execute
.select_cap <- function(liana_res, fun){
    nums <- liana_res %>% map(function(res) nrow(res)) %>% as.numeric
    exec(fun, nums)
}




