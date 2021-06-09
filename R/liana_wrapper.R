#' Run LIANA Wrapper
#'
#' @param seurat_object seurat object
#' @param method method(s) to be run via liana
#' @param resource resource(s) to be used by the methods
#'
#' @import tibble rlang
#' @importFrom purrr map map2
#'
#' @returns A list of method-resource results - i.e. provided resources are run
#' with each method
#'
#' @details LIANA wrapper method that can be used to call each method with
#'  a given set of intercellular resources from the OmniPath universe
liana_wrap <- function(seurat_object, method, resource){
    resource %<>% .select_resource

    .select_method(method) %>%
        map2(names(.),
             function(.method, method_name){
                 if(!(method_name %in% c('squidpy', 'natmi'))){
                     map(resource, function(reso){
                         args <- append(
                             list(seurat_object = seurat_object,
                                  op_resource = reso),
                             options(str_glue('{method_name}.defaults'))[[1]]
                         )
                         exec(.method,  !!!args)
                     })
                 } else{
                     args <- append(
                         list(seurat_object = seurat_object,
                              op_resource = resource),
                         options(str_glue('{method_name}.defaults'))[[1]]
                     )
                     exec(.method,  !!!args)
                 }
             })
}

#' Helper Function to Handle resource choices
#' @param resource names of the resources
.select_resource <- function(resource){
    omni_resources <- readRDS("data/input/omni_resources.rds")

    if(tolower(resource)=="all"){
        omni_resources[as.character(get_lr_resources())]
    } else{
        omni_resources[resource]
    }
}


#' Function to
#' @param method name of the method
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


# Helper Function to  Handle list or not list method calls
.list2tib <- function(res){
    if(length(res)==1){res %>% pluck(1)} else{res}
}
