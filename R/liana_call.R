#' Function to obtain SingleCellSignalR-like scores
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`LRscore`) as calculated
#'    by SingleCellSignalR
get_sca <- function(lr_res,
                    ...){
    liana_call(
        lr_res = lr_res,
        method = "sca",
        ...
    )
}


#' Function to obtain connectome-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`weight_sc`) as calculated
#'    by Connectome
get_connectome <- function(lr_res,
                           ...){
    liana_call(
        lr_res = lr_res,
        method = "connectome",
        ...
        )
}



#' Function to obtain NATMI-like weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with specificity weights (`edge_specificity`)
#'    as calculated by NATMI
get_natmi <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "natmi",
        ...
    )
}



#' Function to obtain logFC weights
#'
#' @inheritParams liana_scores
#' @inheritDotParams liana_call
#'
#' @export
#'
#' @return Returns a tibble with a logFC metric (`logfc_comb`). `logfc_comb` is
#'    calculated as the product of the (1 vs the rest) log2FC for each ligand
#'    and receptor gene
get_logfc <- function(lr_res,
                      ...){
    liana_call(
        lr_res = lr_res,
        method = "logfc",
        ...
    )
}


#' Function to obtain CellPhoneDB-like scores
#'
#' @inheritParams liana_scores
#' @inheritParams cellphonedb_score
#' @inheritDotParams cellphonedb_score
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
get_cellphonedb <- function(lr_res,
                            ...){
    liana_call(
        lr_res = lr_res,
        method = "cellphonedb",
        ...
    )
}


#' Function to obtain CytoTalk-like scores
#'
#' @inheritParams liana_scores
#' @inheritParams cellphonedb_score
#' @inheritDotParams cellphonedb_score
#'
#' @noRd
#'
#' @return Returns a tibble with specificity weights (`pvalue`) as calculated
#'    by CellPhoneDB
get_cytotalk <- function(lr_res,
                         ...){

    liana_call(
        lr_res = lr_res,
        method = "cytotalk",
        ...
    )
}


#' Wrapper Function to obtain scores via liana_pipe
#' @param method name of the method to be called
#' @inheritParams liana_pipe
#' @inheritDotParams liana_pipe
#' @inheritParams liana_scores
#' @inheritParams recomplexify
#'
#' @export
#'
#' @return lr_res modified to be method-specific
liana_call <- function(method,
                       lr_res,
                       ...){

    liana_scores(.score_specs()[[method]],
                 lr_res = lr_res,
                 ...)
}



#' Function to obtain different scoring schemes
#'
#' @param score_object score_object specific to the test obtained from score_specs
#' @param lr_res ligand-receptor DE results and other stats between clusters
#' @param ... dot params passed to `*_score` functions
#'
#' @return lr_res modified to be method-specific
liana_scores <- function(score_object,
                         lr_res,
                         ...){


    # Get expr prop from defaults/kwargs
    expr_prop <- list(...)[["expr_prop"]]
    return_all <- list(...)[["return_all"]]

    supp_columns <- list(...) %>% pluck("supp_columns")
    supp_columns <- union(score_object@add_columns, supp_columns)
    if(return_all){
        supp_columns <- union(supp_columns, c("lrs.to.keep"))
    }

    # join all columns
    all_columns <- c(score_object@columns, supp_columns)
    # Remove NAs for methods that don't have additional columns
    all_columns <- as.character(na.omit(all_columns))

    lr_res %<>%
        select(ligand, receptor,
               ends_with("complex"),
               source, target,
               ends_with("prop"),
               any_of(all_columns))

    lr_res %<>%
        recomplexify(
            lr_res = .,
            columns = score_object@columns,
            add_columns = supp_columns,
            ...)  %>%
        # Select only the relevant columns
        select(source, target,
               ligand.complex, ligand,
               receptor.complex, receptor,
               ends_with("prop"),
               any_of(all_columns)) %>%
        ungroup()

    if(return_all){
        lr_res <- lr_res %>%
            mutate(lrs.to.keep = (receptor.prop >= expr_prop &
                       ligand.prop >= expr_prop))
        rest <- lr_res %>%
            filter(!lrs.to.keep)
        if(score_object@method_name!="cytotalk"){
            lr_res <- lr_res %>%
                filter(lrs.to.keep)
        }
    } else if(!return_all){
        rest <- NULL
    }

    args <-
        append(
            list(lr_res = lr_res,
            score_col = score_object@method_score),
            list(...)
        )

    old_columns <- colnames(lr_res)
    lr_res <- exec(score_object@score_fun, !!!args)
    new_columns <- setdiff(colnames(lr_res), old_columns)

    lr_res %>%
        ungroup() %>%
        .assign_to_filter(lr_res=.,
                          columns = new_columns,
                          expr_prop=expr_prop,
                          return_all = args$return_all,
                          rest=rest) %>%
        # ensure that there are no duplicates (edge cases where multiple subunits
        # have the same expr. - note that we also include method score to ensure
        # that no information is being lost + there are no issues)
        distinct_at(c("source", "target",
                      "ligand.complex", "receptor.complex",
                      score_object@method_score), .keep_all = TRUE)
}


#' Don't filter but assign as worst, and add `lrs.to.keep`
#'
#' @keywords internal
#' @noRd
.assign_to_filter <- function(lr_res,
                              columns,
                              expr_prop,
                              return_all,
                              rest){
    if(!return_all){
        return(lr_res %>%
                   filter(receptor.prop >= expr_prop &
                              ligand.prop >= expr_prop)
               )
    } else{
        lr_res <- bind_rows(lr_res, rest)

        map(columns,function(col){

            # TODO Need to change this
            # Here, I get the newly assigned columns and according to those
            # I set the mins and max...
            # Would require to rework the classes as in Python
            # deal with descending
            if(col %in% c("pvalue", "pval")){
                fun <- "max"
            } else{ # ascending
                fun <- "min"
            }

            set_to <- exec(.fn=fun, lr_res[[col]], na.rm=TRUE)
            lr_res <<- lr_res %>%
                mutate({{ col }} := ifelse(lrs.to.keep, lr_res[[col]], set_to))
        })

        return(lr_res)
    }
}



