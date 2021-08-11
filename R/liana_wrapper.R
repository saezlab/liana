#' Run LIANA Wrapper
#'
#' @param seurat_object seurat object
#' @param method method(s) to be run via liana
#' @param resource resource(s) to be used by the methods (Use `all` to run all resources in one go),
#'   or `custom` to run liana_wrap with an appropriately formatted custom resource, passed via `exernal_resource`
#' @param external_resource external resource in OmniPath tibble format
#' @param .simplify if methods are run with only 1 resource, return a list
#'   of tibbles for each method (default), rather than a list of lists with
#'   method-resource combinations
#' @inheritDotParams liana_defaults
#'
#' @import tibble
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
#'
#' @export
liana_wrap <- function(seurat_object,
                       method = c('natmi', 'connectome', 'logfc',
                                  'cellchat', 'sca', 'squidpy'),
                       resource = c('OmniPath'),
                       external_resource,
                       .simplify = TRUE,
                       ...){
  if(length(setdiff(tolower(method), show_methods())) > 0){
    stop(str_glue("{setdiff(tolower(method), show_methods())} not part of LIANA "))
    }

    if(resource!='custom' & length(setdiff(resource, c(show_resources(), "all"))) > 0){
      stop(str_glue("{setdiff(resource, show_resources())} not part of LIANA "))
      }

  if(resource!='custom'){
    resource %<>% select_resource
  } else{
    resource = list('custom_resource'=external_resource)
  }

  if(any(method %in% c("natmi", "connectome", # change this
                       "logfc", "sca"))){

    lr_results <- resource %>%
      map(function(reso){

        if(is.null(reso)){
          stop("Resource is NULL and LIANA PIPE methods have no default")
        }

        args <- append(
          list("seurat_object" = seurat_object,
               "op_resource" = reso),
          liana_defaults(...)[["liana_pipe"]]
          )

        rlang::invoke(liana_pipe, args)
        }) %>%
      setNames(names(resource))
  }

  .select_method(method) %>%
    map2(names(.),
         safely(function(.method, method_name){
           message(str_glue("Now Running: {stringr::str_to_title(method_name)}"))

           map2(resource, names(resource), function(reso, reso_name){
             # external calls
             if(!(method_name %in% c("natmi", "connectome", # change this
                                     "logfc", "sca"))){
               args <- append(
                 list("seurat_object" = seurat_object,
                      "op_resource" = reso),
                 liana_defaults(...)[[method_name]]
                 )
               rlang::invoke(.method,  args)
               } else {
                 # re-implemented non-permutation approaches
                 args <- append(
                   list("seurat_object" = seurat_object,
                        lr_res = lr_results[[reso_name]]),
                   liana_defaults(...)[["liana_call"]]
                   )
                 rlang::invoke(.method,  args)
                 }
             })
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
        omni_resources[show_resources()]
    } else{
        omni_resources[resource]
    }
}


#' Function to return the appropriate method(s) to be executed
#' @param method name of the method
#'
#' @return A list of method function names (to be called by the LIANA wrapper)
#'
#' @details Adapted from (jvelezmagic); [decoupleR](https://github.com/saezlab/decoupleR/)
#'
#' @noRd
.select_method <- function(method){
    available_method <-
        list(
            # liana_scores
            connectome = expr(get_connectome),
            logfc = expr(get_logfc),
            natmi = expr(get_natmi),
            sca = expr(get_sca),
            # pipes
            squidpy = expr(call_squidpy),
            cellchat = expr(call_cellchat),
            # deprecated
            call_sca = expr(call_sca),
            call_connectome = expr(call_connectome),
            call_natmi = expr(call_natmi),
            call_italk = expr(call_italk)
        )

    method %>%
        tolower() %>%
        match.arg(names(available_method), several.ok = TRUE) %>%
        available_method[.]
}




#' Helper Function to return the methods in LIANA
#'
#' @details methods starting with `call_*` were re-implemented in liana and
#'    albeit their original pipelines (and packages are still supported),
#'    we recommend using the liana re-implementations for efficiency
#'
#' @export
show_methods <- function(){
    c("logfc",
      "natmi",
      "connectome",
      "sca",
      "squidpy",
      "cellchat",
      "call_sca",
      'call_natmi',
      'call_italk',
      'call_connectome')
}

#' Helper Function to return the Resources in LIANA
#'
#' @export
show_resources <- function(){
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
