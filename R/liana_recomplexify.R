#' Helper function to account for complexes in the resources
#'
#' @param lr_cmplx decomplexified* lr_res
#' @param columns columns to account for complexes for
#' @param complex_policy policy how to account for the presence of complexes.
#'   Following the example of \href{https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.ligrec.html}{Squidpy} valid options are:
#'   `'min0'` (default): select the subunit with the change/expression closest to 0
#'    (as in \href{https://github.com/Teichlab/cellphonedb}{CellPhoneDB} when working with expression)
#'   `'all'`: returns all possible subunits separately;
#'   Alternatively, pass any other base or user-defined function that reduces a vector to a single number (e.g. mean, min, etc)
#'
#' @returns complex-accounted lr_res, with both complex and subunit genesymbols
#'
#' @details to be passed before the relevant score_calc function;
#' * decomplexified refers to complexes being broken into subunits which are then
#'  treated as seperate entities and re-assembled into complexes (or 'recomplexified')
#'  by this function
#'
#' @export
#'
#' @importFrom stringr str_split
recomplexify <- function(lr_res,
                         columns,
                         complex_policy){

    if(complex_policy=='all'){ return(lr_res) }

    # Account for Missing Subunits
    recomplexify_env = new.env()
    lr_res %<>% account_missing(recomplexify_env)

    # Account for Complexes
    grps <- c("source", "target",
              "ligand.complex", "receptor.complex")

    columns %>%
        map(function(col){

            col.min <- sym(str_glue("{col}.min"))
            col.flag <- sym(str_glue("{col}.flag"))

            lr_res <<- lr_res %>%
                group_by(across(all_of(grps))) %>%
                mutate( {{ col.min }} := exec(complex_policy, (.data[[col]]))) %>%
                mutate( {{ col.flag }} :=
                            ifelse(.data[[col]]==.data[[col.min]],
                                   TRUE,
                                   FALSE))
        })

    lr_cmplx <- lr_res %>%
        filter(if_all(ends_with("flag"))) %>%
        group_by(across(all_of(c("source", "target",
                                 "ligand", "receptor")))) %>%
        select(source, target,
               ligand.complex, ligand,
               receptor.complex, receptor,
               ends_with("prop"),
               !!columns)

    return(lr_cmplx)
}



#' Function to account for missing subunits
#'
#' @param lr_res LR results as obtained by `lr_pipe`
#' @param env environemnt to which the lr_res will be saved
#'
#' @details this function applies loops over all complexes and assigns 0s to
#'    the stats of any complexes with a missing subunit. The loop is iteratively
#'    applied over the same `lr_res`, to provide a more controlled environment
#'    we pass the environment of recomplexify and modify `lr_res` in place.
account_missing <- function(lr_res, env){

    env$lr_res <- lr_res

    ligand_complexes <- lr_res %>%
        filter(str_detect(.data[["ligand.complex"]], "_")) %>%
        pluck("ligand.complex") %>%
        unique()

    receptor_complexes <- lr_res %>%
        filter(str_detect(.data$receptor.complex, "_")) %>%
        pluck("receptor.complex") %>%
        unique()

    map(receptor_complexes,
        ~missing_subunits_to0(lr_res = lr_res,
                              complex = .x,
                              entity = "receptor",
                              env = env))

    map(receptor_complexes,
        ~missing_subunits_to0(lr_res = lr_res,
                              complex = .x,
                              entity = "ligand",
                              env = env))

    return(env$lr_res)
}


#' Helper Function that assigns 0s to any complexes with missing subunits
#'
#' @param lr_res LR results as obtained by `lr_pipe`
#' @param complex complex of interest
#' @param entity is the complex a 'ligand' or 'receptor'
#'
#' @return A `lr_res` tibble with
missing_subunits_to0 <- function(lr_res, complex, entity, env){

    entity.complex <- str_glue("{entity}.complex")

    complex_split <- # split complex into subunits
        str_split(complex, "_") %>%
        pluck(1)

    # check if subunits are present
    absent_subunits <- setdiff(complex_split,
                               lr_res[[entity]] %>% unique())

    # if there are absent subunits assign 0s and pvalues of 1
    if(length(absent_subunits)>0){
        env$lr_res <- env$lr_res %>%
            # any numeric value to 0
            mutate(across(where(is.numeric) & starts_with(!!entity),
                          ~ ifelse(.data[[entity.complex]]==complex, 0, .))) %>%
            # FDR and pval to 1
            mutate(across(starts_with(!!entity) &
                              (ends_with("FDR") | ends_with("pval")),
                          ~ ifelse(.data[[entity.complex]]==complex, 1, .)))
    }

    return()
}


#' Helper Function which returns the value closest to 0
#' @param vec numeric vector
#'
#' @return value closest to 0
min0 <- function(vec){
    vec[which.min(abs(vec))]
}
