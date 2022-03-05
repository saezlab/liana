#' LIANA wrapper function
#'
#' @param sce `SingleCellExperiment` object or `SeuratObject`
#'
#' @param method method(s) to be run via liana
#'
#' @param resource resource(s) to be used by the methods (`Consensus` by default), Use `all` to run all resources in one go),
#'   or `custom` to run liana_wrap with an appropriately formatted custom resource, passed via `exernal_resource`
#'
#' @param idents the cell identities/labels to be used.
#'
#' @param external_resource external resource in OmniPath tibble format
#'
#' @param .simplify if methods are run with only 1 resource, return a list
#'   of tibbles for each method (default), rather than a list of lists with
#'   method-resource combinations
#'
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
liana_wrap <- function(sce,
                       method = c('natmi', 'connectome', 'logfc',
                                  'sca', 'cellphonedb'),
                       resource = c('Consensus'),
                       idents_col = NULL,
                       external_resource,
                       verbose = TRUE,
                       .simplify = TRUE,
                       ...){

  # Handle object
  sce <- liana_prep(sce, idents_col = idents_col, verbose = verbose)

  # method to lower
  method %<>% stringr::str_to_lower()

  if(resource!='custom' & length(setdiff(resource, c(show_resources(), "all"))) > 0){
    stop(str_glue("{setdiff(resource, show_resources())} not part of LIANA "))
    }

  if(resource!='custom'){
    resource %<>% select_resource # if null OmniPath
  } else{
    resource = list('custom_resource' = external_resource)
  }

  if(any(method %in% c("natmi", "connectome", # change this (logical for internal)
                       "logfc", "sca",
                       "cellphonedb", "cytotalk", "scconnect"))){

    # LIANA pipe map over resource
    lr_results <- resource %>%
      map(function(reso){

        if(is.null(reso)){
          liana_message("Resource was NULL and LIANA's internal methods were run with the `Consensus` resource",
                        verbose = verbose,
                        output = "warning")
          reso <- select_resource("Consensus")[[1]]
        }

        rlang::invoke(liana_pipe,
                      append(
                        list("sce" = sce,
                             "op_resource" =  decomplexify(reso),
                             verbose = verbose),
                        liana_defaults(...)[["liana_pipe"]]
                        )
                      )
        }) %>%
      setNames(names(resource))
  }

  .select_method(method) %>%
    map2(names(.),
         safely(function(.method, method_name){
           liana_message(str_glue("Now Running: {stringr::str_to_title(method_name)}"),
                         verbose = verbose)

           map2(resource, names(resource), function(reso, reso_name){
             if(method_name %in% c("squidpy", "cellchat")){
               # external calls (for complex-informed methods)
               args <- append(
                 list("sce" = sce,
                      "op_resource" = reso),
                 liana_defaults(...)[[method_name]]
                 )
               rlang::invoke(.method, args)

             } else if(startsWith(method_name, "call_")){
               # external calls (for methods who don't deal with complexes)
               args <- append(
                 list(sce = sce,
                      "op_resource" = reso %>% {if (!is.null(reso)) decomplexify(.) else .}),
                 liana_defaults(...)[[method_name]]
               )
               rlang::invoke(.method, args)

             } else if(method_name == "cellphonedb"){
               # get lr_res for this specific resource
               lr_res <- .filt_liana_pipe(lr_results[[reso_name]],
                                          method_name,
                                          ...)

               # permutation-based approaches
               perm_means <-
                 rlang::invoke(
                   get_permutations,
                   append(
                     list(lr_res = lr_res,
                          sce = sce,
                          verbose = verbose
                          ),
                     liana_defaults(...)[["permutation"]]
                     ))

               args <- pmap( # append multiple lists
                 list(
                   list(
                     list(lr_res = lr_res,
                          perm_means = perm_means,
                          verbose = verbose),
                     liana_defaults(...)[[method_name]],
                     liana_defaults(...)[["liana_call"]]
                     )
                   ), c) %>%
                 flatten

               rlang::invoke(.method, args)

             } else if(method_name == "cytotalk"){
               lr_res <- .filt_liana_pipe(lr_results[[reso_name]],
                                          method_name,
                                          ...)

               args <- pmap( # append multiple lists
                 list(
                   list(
                     list(
                       lr_res = lr_res,
                       sce = sce
                     ),
                     liana_defaults(...)[[method_name]],
                     liana_defaults(...)[["liana_call"]]
                     )
                 ), c) %>%
                 flatten

               rlang::invoke(.method, args)

            } else {
              # re-implemented non-permutation approaches
              args <- append(
                list(
                  "sce" = sce,
                  lr_res = .filt_liana_pipe(lr_results[[reso_name]],
                                            method_name,
                                            ...)
                  ),
                liana_defaults(...)[["liana_call"]]
                )
              rlang::invoke(.method, args)
              }
             })
           }, quiet = FALSE)) %>%
    # format errors
    {
      `if`(.simplify,
           map(., function(elem)
             .list2tib(.list2tib(compact(elem)))),
           .)
    } %>% .list2tib()
}


#' Helper Function to Handle resource choices
#'
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
            scconnect = expr(get_scconnect),
            # liana_permutes
            cellphonedb = expr(get_cellphonedb),
            # cytotalk
            cytotalk = expr(get_cytotalk),
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
      'call_connectome',
      "cellphonedb"
      )
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


#' Helper Function to do prop filtering of liana_pipe output
#' @param liana_res resutls from `liana_pipe`
#' @param method current method
#' @inheritDotParams liana_defaults
.filt_liana_pipe <- function(liana_res,
                             method,
                             ...){

  # Convert prop_filt logical to 0 or 1
  #  (step not required, but should be explicit)
  filt_or_not <- as.integer(liana_defaults(...)[[method]]$prop_filt)

  # set to value of expr_prop from defaults, or 0 if not to filter
  expr_prop <- liana_defaults(...)$expr_prop * filt_or_not

  if(expr_prop > 0){
    liana_res %>%
      filter(receptor.prop >= expr_prop & ligand.prop >= expr_prop)
  } else if(expr_prop == 0){
    liana_res
  } else{ # shouldn't happen (sanity check)
    stop("Expression Proportion is Negative")
  }

}

