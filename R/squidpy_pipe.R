#' Call Squidpy Pipeline via reticulate with OmniPath and format results
#' @param seurat_object Seurat object as input
#' @param omni_resources List of OmniPath resources
#' @param python_path path to python version to use in reticulate
#' @param .seed used to python seed
#' @returns A list of Squidpy results for each resource
#' @details CellPhoneDB v2 algorithm implementation in Python
#' Stats:
#' Mean expr
#' pval from shuffled cluster
#' @import reticulate tibble
#' @export
call_squidpyR <- function(seurat_object,
                          omni_resources,
                          python_path,
                          .seed = 1004,
                          .ident = "seurat_annotations"){

    reticulate::use_python(python_path)
    py$pd <- reticulate::import("pandas")

    if("DEFAULT" %in% toupper(names(omni_resources))){ # to be replaced
        # omni_resources$Default <- omni_resources$CellPhoneDB
        omni_resources$Default <- NULL
    }

    op_resources <- map(omni_resources, function(x) x %>%
                              select(
                                  uniprot_source = source,
                                  unprot_target = target,
                                  source = source_genesymbol,
                                  target = target_genesymbol,
                                  category_intercell_source,
                                  category_intercell_target
                                  )) %>%
        unname() # unname list, so passed as list not dict to Python

    # Call Squidpy
    reticulate::source_python("R/squidpy_pipe.py")
    py_set_seed(.seed)

    py$squidpy_results <- py$call_squidpy(op_resources, # to include **kwargs
                                          GetAssayData(seurat_object), #expr
                                          seurat_object[[]], # meta
                                          GetAssay(seurat_object)[[]], # feature_meta
                                          .ident)

    squidpy_pvalues <- py$squidpy_results$pvalues %>% setNames(names(omni_resources))
    squidpy_means <- py$squidpy_results$means %>% setNames(names(omni_resources))
    squidpy_metadata <- py$squidpy_results$meta %>% setNames(names(omni_resources))


    squidpy_results <- map(names(omni_resources),
                           function(x)
                               squidpy_reformat(.name=x,
                                                .pval_list = squidpy_pvalues,
                                                .mean_list = squidpy_means,
                                                .meta_list = squidpy_metadata)) %>%
        setNames(names(omni_resources)) %>%
        map(function(res) res %>%
                select(1:3, means, pvalue, uniprot_source, unprot_target) %>%
                rename(ligand = source,
                       receptor = target) %>%
                separate(pair, sep = "_", into=c("source", "target")))

    return(squidpy_results)
}


#' Helper function to reformat squidpy function r
#' @param .name omnipath resource name
#' @param .pval_list p-value results from different dbs as a list from squidpy
#' @param .mean_list mean list from squidpy
#' @export
squidpy_reformat <- function(.name,
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
