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
                            cap = NULL,
                            get_ranks = TRUE){

    if(!is_tibble(liana_res[[1]]) && is.null(resource)){
        stop("Please provide provide a name for the resource, ",
                 "otherwise the first resource will be plucked!")
    } else if(!is_tibble(liana_res[[1]])){
        liana_res %<>% map(function(m_results) m_results %>% pluck(resource))
    }

    cap %<>% `%||%`(.select_cap(liana_res, set_cap))

    liana_res %>%
        map2(names(.), function(res, method_name){
            method_score <- .rank_specs()[[method_name]]@method_score
            desc_order <- .rank_specs()[[method_name]]@descending_order

            .method = sym(as.character(str_glue("{method_name}.{method_score}")))
            .rank_col = sym(as.character(str_glue("{method_name}.rank")))

            res %>%
                top_n(n=if_else(desc_order,
                                cap,
                                -cap),
                      wt=!!sym(method_score)) %>%
                mutate( {{ .rank_col }} := .rank_enh(.data[[method_score]], desc_order)) %>%
                arrange(!!.rank_col) %>%
                rename( {{ .method }} := method_score) %>%
                select(source, ligand, target, receptor, !!.method, !!.rank_col) %>%
                distinct() %>%
                as_tibble()
        }) %>%
        purrr::reduce(., full_join, by = c("source", "ligand", # Join all res
                                           "target", "receptor")) %>%
        {`if`(get_ranks, .liana_consensus(. ,cap))}
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
#' @param liana_agg Aggregated method results
#' @param cap Value assigned to NA
#'
#' @return aggregated liana tibble with consensus ranks
.liana_consensus <- function(liana_agg, cap){
    liana_agg %>%
        mutate_at(vars(ends_with(".rank")),
              ~ replace(., is.na(.), cap)) %>% # assign .rank_cap to NA
        mutate(mean_rank = pmap_dbl(select(., ends_with(".rank")),
                                    function(...) mean(c(...))),
               median_rank = pmap_dbl(select(., ends_with(".rank")),
                                      function(...) median(c(...))))  %>%
        arrange(mean_rank) %>%
        select(source, ligand, target, receptor,
               ends_with("_rank"), everything())
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



#' S4 Class used to generate aggregate/consesus scores for the methods.
#'
#' @name RankSpecifics-class
#'
#' @field method_name name of the method (e.g. cellchat)
#' @field method_score The interaction score provided by the method (typically
#' the score that reflects the specificity of interaction)
#' @field descending_order whether the score should be interpreted in
#'  descending order (i.e. highest score for an interaction is most likely)
#'
#' @exportClass RankSpecifics
setClass("RankSpecifics",
         slots=list(method_name="character",
                    method_score="character",
                    descending_order="logical"))

#' Rank Specs Holder
#'
#' @return list of RankSpecifics objects for each method
#'
#' @noRd
.rank_specs <- function(){
    list(
        "cellchat" =
            methods::new(
                "RankSpecifics",
                method_name = "cellchat",
                method_score = "pval",
                descending_order = FALSE
            ),
        "connectome" =
            methods::new(
                "RankSpecifics",
                method_name = "connectome",
                method_score = "weight_sc",
                descending_order = TRUE
            ),
        "italk" =
            methods::new(
                "RankSpecifics",
                method_name = "italk",
                method_score = "weight_comb",
                descending_order = TRUE
            ),
        "natmi" =
            methods::new(
                "RankSpecifics",
                method_name = "natmi",
                method_score = "edge_specificity",
                descending_order = TRUE
            ),
        "sca" = methods::new(
            "RankSpecifics",
            method_name = "sca",
            method_score = "LRscore",
            descending_order = TRUE
        ),
        "squidpy" =
            methods::new(
                "RankSpecifics",
                method_name = "Squidpy",
                method_score = "pvalue",
                descending_order = FALSE
            )
    )
}

