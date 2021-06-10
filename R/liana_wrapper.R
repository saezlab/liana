#' Run LIANA Wrapper
#'
#' @param seurat_object seurat object
#' @param method method(s) to be run via liana
#' @param resource resource(s) to be used by the methods
#' @param .simplify if methods are run with only 1 resource, return a list
#'   of tibbles for each method, rather than a list of lists with
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
#'
#' @details LIANA wrapper method that can be used to call each method with
#'  a given set of intercellular resources from the OmniPath universe
#'
#' @export
liana_wrap <- function(seurat_object,
                       method,
                       resource,
                       .simplify = TRUE,
                       ...){
    resource %<>% .select_resource

    .select_method(method) %>%
        map2(names(.),
             safely(function(.method, method_name){
                 if(!(method_name %in% c('squidpy', 'natmi'))){
                     map(resource, function(reso){
                         args <- append(
                             list(seurat_object = seurat_object,
                                  op_resource = reso),
                             liana_defaults(...)[[method_name]]
                         )
                         exec(.method,  !!!args)
                     }) %>% {`if`(.simplify, .list2tib(.)) }
                 } else{
                     args <- append(
                         list(seurat_object = seurat_object,
                              op_resource = resource),
                         liana_defaults(...)[[method_name]]
                     )
                     exec(.method,  !!!args) %>% {`if`(.simplify, .list2tib(.)) }
                 }
             }, quiet = FALSE)) %>%
        map(function(elem) .list2tib(compact(elem))) # format result/errors
}

#' Helper Function to Handle resource choices
#' @param resource names of the resources
.select_resource <- function(resource){
    omni_resources <-
        readRDS(system.file(package = 'liana', "omni_resources.rds"))

    if(tolower(resource)=="all"){
        omni_resources[as.character(get_lr_resources())]
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
.select_method <- function(method){
    available_method <-
        list(
            cellchat = expr(call_cellchat),
            connectome = expr(call_connectome),
            italk = expr(call_italk),
            natmi = expr(call_natmi),
            sca = expr(call_sca),
            squidpy = expr(call_squidpyR)
        )

    method %>%
        tolower() %>%
        match.arg(names(available_method), several.ok = TRUE) %>%
        available_method[.]
}


#' Helper Function to  Handle list or not list method calls
#' @param res list of lists (e.g. Method-Resource Results from the LIANA pipe)
#'
#' @return The first element of the list
.list2tib <- function(res){
    if(length(res)==1){res %>% pluck(1)} else{res}
}
