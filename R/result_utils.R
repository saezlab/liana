#' Get top hits
#' @param spec_list list of spec objects with ligrec results
#' @return A list of top hits per tool/tool_parameter
get_top_hits <- function(spec_list, .tn=200){
    names(spec_list) %>%
        map(function(method_name){

            parnams <- names(spec_list[[method_name]]@method_scores)

            map(parnams, function(parm){
                spec_list[[method_name]]@method_results %>%
                    map(function(res){
                        res %>%
                            distinct() %>%
                            top_n(n=if_else(spec_list[[method_name]]@method_scores[[parm]],
                                            .tn,
                                            -.tn),
                                  wt=!!sym(parm)) %>%
                            as_tibble()
                    })
            }) %>%
                {if(length(spec_list[[method_name]]@method_scores) > 1)
                    setNames(., str_glue("{method_name}_{parnams}"))
                    else setNames(., method_name)
                }
        }) %>% setNames(names(spec_list)) %>%
        flatten() %>%
        setNames(str_replace_all(names(.),
                                 pattern = "_",
                                 replacement = "\\."))
}


#' S4 Class used to format benchmark output.
#' @name MethodSpecifics-class
#'
#' @field method_name
#' @field method_results
#' @field method_scores
#'
#' @exportClass MethodSpecifics
setClass("MethodSpecifics",
         slots=list(method_name="character",
                    method_results = "list",
                    method_scores="list"))
