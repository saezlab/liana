#' Call Squidpy Pipeline via reticulate with OmniPath and format results
#' @param seurat_object Seurat object as input
#' @param op_resource List of OmniPath resources
#' @param .seed used to python seed
#' @param conda_env python conda environment to run Squidpy; set to liana_env by default
#' @param ... kwargs passed to Squidpy; For more information see:
#'   \link{https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.ligrec.html#squidpy.gr.ligrec}
#'
#' @import reticulate tibble
#' @importFrom tidyr pivot_longer
#' @importFrom Seurat GetAssay GetAssayData
#' @importFrom dplyr left_join
#'
#' @details CellPhoneDB v2 algorithm implementation in Python
#'
#' @returns A list of Squidpy results for each resource
#'
#' @export
call_squidpy <- function(seurat_object,
                          op_resource,
                          .seed = 1004,
                          conda_env = NULL,
                          ...){

    kwargs <- list(...)

    reticulate::use_condaenv(condaenv = conda_env %>% `%||%`("liana_env"),
                             conda = "auto",
                             required = TRUE)
    py$pd <- reticulate::import("pandas")

    if("DEFAULT" %in% toupper(names(op_resource))){
        op_resource$Default <- NULL

    }

    op_resources <- map(op_resource, function(x) x %>%
                              select(
                                  uniprot_source = source,
                                  unprot_target = target,
                                  source = source_genesymbol,
                                  target = target_genesymbol,
                                  category_intercell_source,
                                  category_intercell_target
                                  )
                        ) %>%
        unname() # unname r list, so its passed as list to Python

    # Call Squidpy
    reticulate::source_python(system.file(package = 'liana', "squidpy_pipe.py"))
    py_set_seed(.seed)

    py$squidpy_results <- py$call_squidpy(op_resources,
                                          GetAssayData(seurat_object), #expr
                                          seurat_object[[]], # meta
                                          GetAssay(seurat_object)[[]], # feature_meta
                                          kwargs # passed to squidpy.gr.ligrec
                                          )

    squidpy_pvalues <- py$squidpy_results$pvalues %>% setNames(names(op_resource))
    squidpy_means <- py$squidpy_results$means %>% setNames(names(op_resource))
    squidpy_metadata <- py$squidpy_results$meta %>% setNames(names(op_resource))


    squidpy_results <- map(names(op_resource),
                           function(x)
                               FormatSquidpy(.name=x,
                                                .pval_list = squidpy_pvalues,
                                                .mean_list = squidpy_means,
                                                .meta_list = squidpy_metadata)) %>%
        setNames(names(op_resource)) %>%
        map(function(res) res %>%
                select(1:3, means, pvalue, uniprot_source, unprot_target) %>%
                rename(ligand = source,
                       receptor = target) %>%
                separate(pair, sep = "_", into=c("source", "target"))) %>%
        .list2tib()

    return(squidpy_results)
}


#' Helper function to reformat Squidpy function r
#' @param .name omnipath resource name
#' @param .pval_list p-value results from different dbs as a list from squidpy
#' @param .mean_list mean list from squidpy
#'
#' @importFrom dplyr left_join
#'
#' @noRd
FormatSquidpy <- function(.name,
                             .pval_list,
                             .mean_list,
                             .meta_list){
    x_pval <- .pval_list[[.name]] %>%
        py_to_r() %>%
        pivot_longer(cols = 3:ncol(.),
                     values_to="pvalue",
                     names_to="pair")

    x_mean <- .mean_list[[.name]] %>%
        py_to_r() %>%
        pivot_longer(cols = 3:ncol(.),
                     values_to = "means",
                     names_to = "pair")

    x_meta <- .meta_list[[.name]] %>% py_to_r()

    res_formatted <- left_join(x_pval, x_mean,
                               by = c("source", "target", "pair")) %>%
        left_join(x_meta, by = c("source", "target"))

    return(res_formatted)
}
