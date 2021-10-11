#' Call Squidpy Pipeline via reticulate with OmniPath and format results
#' @param seurat_object Seurat object as input
#' @param op_resource Tibble or list of OmniPath resources, typically obtained via
#'    \code{\link{select_resource}}
#' @param seed seed passed to squidpy's ligrec function
#' @param conda_env python conda environment to run Squidpy; set to liana_env by default
#' @param ... kwargs passed to Squidpy; For more information see:
#'   \link{https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.ligrec.html#squidpy.gr.ligrec}
#' @param expr_prop minimum proportion of gene expression per cell type (0 by default),
#'  yet perhaps one should consider setting this to an appropriate value between 0 and 1,
#'  as an assumptions of these method is that communication is coordinated at the cluster level.
#'
#' @import reticulate tibble
#' @importFrom tidyr pivot_longer
#' @importFrom Seurat GetAssay GetAssayData
#' @importFrom dplyr left_join
#'
#' @details CellPhoneDBv2 algorithm re-implementation in Python.
#'
#' Note that `cluster_key` is a parameter passed to the Squidpy function,
#' by default this will be set to the default Ident of the Seurat object.
#'
#' @returns A list of Squidpy results for each resource
#'
#' @export
call_squidpy <- function(seurat_object,
                         op_resource,
                         seed = 1004,
                         conda_env = NULL,
                         assay.type = "logcounts",
                         assay = "RNA",
                         ...){

    # required until I make the interchaengeable Seurat/SCE
    if(assay.type=="logcounts"){
        assay.type = "data"
    }

    if(is_tibble(op_resource)){
        op_resource <- list("placeholder" = op_resource)
    }

    kwargs <- list(...)
    kwargs$cluster_key %<>% `%||%`(.get_ident(seurat_object))
    kwargs$seed <- as.integer(seed)

    if(length(kwargs$cluster_key) == 0){
        stop("Squidpy: Cluster annotations missing! Please specificy a column")
    } else{
        message(str_glue("Squidpy: running with {kwargs$cluster_key} as cluster_key"))
    }

    if(!is.factor(seurat_object@meta.data[[kwargs$cluster_key]])){
        seurat_object@meta.data[[kwargs$cluster_key]] <-
            seurat_object@meta.data %>%
            pull(kwargs$cluster_key) %>%
            as.factor()
        message(str_glue("Squidpy: {kwargs$cluster_key} was converted to factor"))
    }

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
    py_set_seed(seed)

    # Check if assay meta.features match object rownames
    # if not assign a placeholder (Seurat-specific fix)
    if(nrow(GetAssay(seurat_object)[[]]) != nrow(seurat_object) |
       ncol(seurat_object@assays[[assay]]@meta.features)==0){
        message("Meta features were reassigned to match object")

        seurat_object@assays[[assay]]@meta.features <-
            data.frame(row.names = rownames(seurat_object),
                       placeholder = rownames(seurat_object))
    }

    py$squidpy_results <- py$call_squidpy(op_resources,
                                          GetAssayData(seurat_object,
                                                       slot=assay.type), # raw expr
                                          seurat_object[[]], # meta
                                          GetAssay(seurat_object,
                                                   assay=assay)[[]], # feature_meta
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
                rename(ligand = source,
                       receptor = target) %>%
                separate(pair, sep = "_", into=c("source", "target")) %>%
                select(source, target,
                       ligand, receptor,
                       means, pvalue)) %>%
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


#' Helper Function to get active Ident (cluster annotation column) from the
#'   Seurat object metadata
#'
#' @param seurat_object seurat object with metadata information
#'
#' @noRd
.get_ident <- function(seurat_object){
    map(names(seurat_object@meta.data),
        function(x){
            p <- seurat_object@meta.data %>%
                select(sym(x)) %>%
                rownames_to_column("names") %>%
                deframe()

            if(identical(p, Seurat::Idents(seurat_object))){
                return(x)
            }

            return()
        }) %>% compact %>%
        as.character %>%
        pluck(1) # to handle scenario when there are two identical columns
}
