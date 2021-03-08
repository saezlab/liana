#' Call Squidpy Pipeline via reticulate with OmniPath and format results
#' @param seurat_object Seurat object as input
#' @param omni_resources List of OmniPath resources
#' @param python_path path to python version to use in reticulate
#' @param .seed used to python seed
#' @returns A list of Squidpy results for each resource
#' @details CellPhoneDB v2 algorithm implementation in Python
#' Stats:
#' Mean expr
#' pval
call_squidpyR <- function(seurat_object,
                          omni_resources,
                          python_path,
                          .seed = 1004,
                          ident = "seurat_annotations"){

    # prep seurat data for transfer
    exprs <- GetAssayData(seurat_object)
    meta <- seurat_object[[]]
    feature_meta <- GetAssay(seurat_object)[[]]
    embedding <- Embeddings(seurat_object, "umap")


    reticulate::use_python(python_path)
    py$pd <- reticulate::import("pandas")

    # Call Squidpy
    reticulate::source_python("scripts/pipes/squidpy_pipe.py")
    py_set_seed(.seed)
    py$squidpy_results <- py$call_squidpy(names(omni_resources),
                                          exprs,
                                          meta,
                                          feature_meta,
                                          embedding,
                                          ident)

    squidpy_pvalues <- py$squidpy_results$pvalues %>% setNames(names(omni_resources))
    squidpy_means <- py$squidpy_results$means %>% setNames(names(omni_resources))


    squidpy_results <- map(names(omni_resources),
                           function(x)
                               squidpy_reformat(.name=x,
                                                .pval_list = squidpy_pvalues,
                                                .mean_list = squidpy_means)) %>%
        setNames(names(omni_resources)) %>% #*
        # swap positions for means and pvalue
        map(function(x) x %>% select(1:3, means, pvalue) %>%
                rename(ligand = source,
                       receptor = target) %>%
                separate(pair, sep = "_", into=c("source", "target")))


    return(squidpy_results)
}


#' Helper function to reformat squidpy function r
#' @param .name omnipath resource name
#' @param .pval_list p-value results from different dbs as a list from squidpy
#' @param .mean_list mean list from squidpy
squidpy_reformat <- function(.name,
                             .pval_list,
                             .mean_list){
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

    res_formatted <- left_join(x_pval, x_mean)
    return(res_formatted)
}

