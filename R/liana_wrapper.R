#' Run LIANA Wrapper
#'
#' @param seurat_object seurat object
#' @param method method(s) to be run via liana
#' @param resource resource(s) to be used by the methods
#' @param .simplify if methods are run with only 1 resource, return a list
#'   of tibbles for each method (default), rather than a list of lists with
#'   method-resource combinations
#' @param ... Pass custom method parameters to the methods via the liana wrapper
#'    See \link{liana::liana_defaults} for more information
#'
#' @import tibble rlang
#' @importFrom magrittr %<>% %>%
#' @importFrom purrr map map2 safely compact
#'
#' @returns A list of method-resource results - i.e. provided resources are run
#' with each method
#' If only one resource is selected, a single tibble (with results for that
#'  resource) will be returned for each of the selected methods
#'
#' @details LIANA wrapper method that can be used to call each method with
#'  a given set of intercellular resources from the OmniPath universe.
#' Use `all` to run all resources in one go.
#'
#' @export
liana_wrap <- function(seurat_object,
                       method = c('cellchat', 'connectome', 'italk',
                                  'natmi', 'sca', 'squidpy'),
                       resource = c('OmniPath'),
                       .simplify = TRUE,
                       ...){

    if(length(setdiff(tolower(method), .get_methods())) > 0){
        stop(str_glue("{setdiff(method, .get_methods())} not part of LIANA "))
    }

    if(length(setdiff(resource, c(get_resources(), "all"))) > 0){
        stop(str_glue("{setdiff(resource, get_resources())} not part of LIANA "))
    }

    resource %<>% select_resource

    .select_method(method) %>%
        map2(names(.),
             safely(function(.method, method_name){
                 if(!(method_name %in% c('squidpy', 'natmi'))){
                     map(resource, function(reso){
                         args <- append(
                             list("seurat_object" = seurat_object,
                                  "op_resource" = reso),
                             liana_defaults(...)[[method_name]]
                         )
                         exec(.method,  !!!args)
                     })
                 } else{
                     args <- append(
                         list("seurat_object" = seurat_object,
                              "op_resource" = resource),
                         liana_defaults(...)[[method_name]]
                     )
                     exec(.method,  !!!args)
                 }
             }, quiet = FALSE)) %>%
        # format errors
        {`if`(.simplify, map(., function(elem)
            .list2tib(.list2tib(compact(elem)))))}

}

#' Helper Function to Handle resource choices
#' @param resource names of the resources
#'
#' @details This function simply reads omni_resources.rds and returns the resources.
#'    Any of the resources can also be obtained via the same file.
#'    or the `compile_ligrec` function, which querries and assmelbes the
#'    resources via `OmniPathR`.
#'
#'    `Default` - The Default (inbuilt) resource for each of the methods;
#'    if using the `call_*` functions, the default resource is used by
#'    passing *NULL* to the resource parameter.
#'    `Reshuffled` - a reshuffled (randomized control) version of ConnectomeDB
#'
#' @export
select_resource <- function(resource){
    omni_resources <-
        readRDS(system.file(package = 'liana', "omni_resources.rds"))

    if(tolower(resource)=="all"){
        omni_resources[get_resources()]
    } else{
        omni_resources[resource]
    }
}


#' Function to return the appropriate method(s) to be executed
#' @param method name of the method
#'
#' @return A list of method function names (to be called by the LIANA wrapper)
#'
#' @details Adapted from (jvelezmagic); decoupleR\https://github.com/saezlab/decoupleR/
#'
#' @noRd
.select_method <- function(method){
    available_method <-
        list(
            cellchat = expr(call_cellchat),
            connectome = expr(call_connectome),
            italk = expr(call_italk),
            natmi = expr(call_natmi),
            sca = expr(call_sca),
            squidpy = expr(call_squidpy)
        )

    method %>%
        tolower() %>%
        match.arg(names(available_method), several.ok = TRUE) %>%
        available_method[.]
}




#' Helper Function to return the methods in LIANA
#' @noRd
.get_methods <- function(){
    c("italk", "squidpy", "natmi", "cellchat", "connectome", "sca")
}

#' Helper Function to return the Resources in LIANA
#'
#' @export
get_resources <- function(){
    as.character(names(
        readRDS(system.file(package = 'liana', "omni_resources.rds"))
        ))
}


#' Helper Function to  Handle list or not list method calls
#' @param res list of lists (e.g. Method-Resource Results from the LIANA pipe)
#' @param .x name of the element to pluck
#'
#' @return The first element of the list
#'
#' @noRd
.list2tib <- function(res, .x = NULL){
    if(length(res)==1){
        .x %<>% `%||%`(1)
        res %>% pluck(.x)
    } else if(is.character(.x)){
        res %>% pluck(.x)
    } else{ res }
}
