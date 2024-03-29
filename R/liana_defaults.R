#' Function to pass Default Arguments for each method
#'
#'
#' @param expr_prop minimum proportion of gene expression per cell type (0.1 by default).
#' Note that when working with complexes, the minimum subunit proportion will
#' be used for filtering.
#' @param complex_policy policy how to account for the presence of complexes.
#'
#' @param seed random seed integer
#' @param parallelize whether to parallelize cellphonedb-like
#' @param workers number of workers to be called
#'
#' @param liana_pipe.params list of Parameters passed to NATMI \code{\link{liana_pipe}}
#' @param liana_call.params list of Parameters passed to NATMI \code{\link{liana_call}}
#' @param cellchat.params list of Parameters passed to CellChat \code{\link{call_cellchat}}
#' @param squidpy.params list of Parameters passed to Squidpy \code{\link{call_squidpy}}
#' @param call_connectome.params list of Parameters passed to Connectome \code{\link{call_connectome}}
#' @param call_italk.params list of Parameters passed to iTALK \code{\link{call_italk}}
#' @param call_natmi.params list of Parameters passed to NATMI \code{\link{call_natmi}}
#' @param call_sca.params list of Parameters passed to SingleCellSignalR \code{\link{call_sca}}
#' @param assay Assay name passed to `call_italk`, `call_sca`, `call_cellchat`,
#'    and `call_connectome`
#' @param cellphonedb.params list of Parameters passed to liana's internal
#' cellphonedb implementation \code{\link{cellphonedb_score}}
#' @param natmi.params list of Parameters passed to liana's internal
#' edge_specificity implementation \code{\link{natmi_score}}
#' @param sca.params list of Parameters passed to liana's internal
#' LRScore implementation \code{\link{sca_score}}
#' @param connectome.params list of Parameters passed to liana's internal
#' connectome's weight_sc implementation \code{\link{connectome_score}}
#' @param cytotalk.params list of Parameters passed to liana's internal
#' crosstalk scores implementation \code{\link{cytotalk_score}}
#' @param logfc.params list of Parameters passed to liana's internal
#' logFC implementation \code{\link{logfc_score}}
#' @param permutation.params list of parameters passed to permutation methods
#' @inheritParams liana_pipe
#'
#' @details The default parameters for each method can also be overwritten by
#'  manually passing a list of parameters for the appropraite method
#'   \code{\link{liana_wrap}}
#'
#' Further, each `get_*` method will by default obtain the default params passed
#'    via \code{\link{liana_pipe}} and \code{\link{liana_call}}. This is done so that most steps
#'    required for the calculation of these methods are undertaken only once.
#'
#' NB! LIANA's internal methods are made consistent. There is no reason to pass
#' specific parameters to any of them. Thus, it is best that one sticks to the
#' non-nested parameters of this function (i.e. excluding `.params`),
#' unless a very specific reason requires any of LIANA's internal parameters
#' to be changed.
#'
#'
#' @return A list of the default parameters for each method
#'
#' @export
#'
#' @examples
#' liana_path <- system.file(package = "liana")
#' # load testdata
#' testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
#' # get a `named` list with all default parameters passed to liana.
#' def_params <- liana_defaults()
#' # any of these can then be overwritten and passed to `...` in `liana_wrap`
#' # with the `.params` suffix to the parameter name type. For example,
#' liana_res <- liana_wrap(testdata,
#'                         permutation.params = list(nperms=2),
#'                         liana_pipe.params = list(test.type='wilcox'))
liana_defaults <- function(
    assay = "RNA",
    assay.type = "logcounts",
    expr_prop = 0.1,
    seed = 1004,
    complex_policy = "mean0",
    parallelize = FALSE,
    workers = 8,
    permutation.params = NULL,
    liana_pipe.params = NULL,
    liana_call.params = NULL,
    cellphonedb.params = NULL,
    natmi.params = NULL,
    sca.params = NULL,
    connectome.params = NULL,
    cytotalk.params = NULL,
    logfc.params = NULL,
    cellchat.params = NULL,
    squidpy.params = NULL,
    call_sca.params = NULL,
    call_natmi.params = NULL,
    call_connectome.params = NULL,
    call_italk.params = NULL,
    ...){

    # Internal ----
    # LIANA_pipe defaults
    liana_pipe.defaults <- list(
        test.type = "wilcox",
        pval.type = "all",
        assay = assay,
        assay.type = assay.type
    )

    # liana_call.defaults
    liana_call.defaults <- list(
        complex_policy = complex_policy,
        expr_prop = expr_prop
    )

    # Permutation defaults (permutations to be used in e.g. CPDB)
    permutation.defaults <- list(
        nperms = 1000,
        parallelize = parallelize,
        workers = workers,
        seed = seed
    )

    # Define Defaults
    # CellPhoneDB defaults
    cellphonedb.defaults <- list(
        workers = workers,
        parallelize = parallelize,
        prop_filt = TRUE
        )

    # NATMI
    natmi.defaults <- list(
        prop_filt = TRUE
    )

    # logFC
    logfc.defaults <- list(
        prop_filt = TRUE
    )


    ## Connectome and Cytotalk calculate scores are calculated at
    ## the cell-cluster-pair level -> We don't apply prop filtering to those
    # Connectome
    connectome.defaults <- list(
        prop_filt = TRUE
    )

    # CytoTalk
    cytotalk.defaults <- list(
        assay.type = assay.type,
        seed = seed,
        prop_filt = FALSE
    )

    # SCA
    sca.defaults <- list(
        prop_filt = TRUE
    )


    # External ----

    # Squidpy defaults
    squidpy.defaults <- list(
        cluster_key = NULL,
        n_perms = 1000,
        threshold = expr_prop,
        seed = as.integer(seed),
        assay = assay,
        assay.type = assay.type
    )

    # CellChat Default
    cellchat.defaults <- list(
        nboot = 100,
        expr_prop = expr_prop,
        exclude_anns = NULL,
        thresh = 1,
        assay = assay,
        .normalize = FALSE,
        .do_parallel = FALSE,
        .raw_use = TRUE,
        organism = "human",
        de_thresh = 0.05
    )

    # SingleCellSignalR Defaults
    call_sca.defaults <- list(
        assay = assay,
        .format = TRUE,
        s.score = 0,
        logFC = log2(1.5)
        )

    # Connectome Defaults
    call_conn.defaults <- list(
        min.cells.per.ident = 1,
        p.values = TRUE,
        calculate.DOR = FALSE,
        assay = assay,
        .format = TRUE,
        .spatial = FALSE
        )

    # NATMI Defaults
    call_natmi.defaults <- list(
        expr_file = "test_em.csv",
        meta_file = "metadata.csv",
        output_dir = "NATMI_results",
        assay = assay,
        assay.type = "logcounts",
        num_cor = 4,
        .format = TRUE,
        .overwrite_data = TRUE,
        .seed = seed,
        .natmi_path = NULL,
        .delete_input_output = FALSE,
        reso_name = "placeholder"
        )

    # iTalk Defaults
    call_italk.defaults <- list(
        assay = assay,
        .format = TRUE,
        .DE = TRUE
        )

    # List of Defaults (re-assigned if needed)
    default_args <- list(
        "expr_prop" = expr_prop,

        # liana_scores (passed to get_* functions)
        "permutation" = permutation.params %<>%
            reassign_params(., permutation.defaults),

        "liana_pipe" = liana_pipe.params %<>%
            reassign_params(., liana_pipe.defaults),

        # this thing needs to be either completely removed or moved to liana_wrap
        "liana_call" = liana_call.params %<>%
            reassign_params(., liana_call.defaults),

        "cellphonedb" = cellphonedb.params %<>%
            reassign_params(., cellphonedb.defaults),

        "cytotalk" = cytotalk.params %<>%
            reassign_params(., cytotalk.defaults),

        "connectome" = connectome.params %<>%
            reassign_params(., connectome.defaults),

        "natmi" = natmi.params %<>%
            reassign_params(., natmi.defaults),

        "sca" = sca.params %<>%
            reassign_params(., sca.defaults),

        "logfc" = logfc.params %<>%
            reassign_params(., logfc.defaults),

        # external methods
        "call_cellchat" = cellchat.params %<>%
            reassign_params(., cellchat.defaults),

        'call_squidpy' = squidpy.params %<>%
            reassign_params(., squidpy.defaults),

        # call_* functions/pipes
        'call_sca' = call_sca.params %<>%
            reassign_params(., call_sca.defaults),

        'call_connectome' = call_connectome.params %<>%
            reassign_params(., call_conn.defaults),

        'call_natmi' = call_natmi.params %<>%
            reassign_params(., call_natmi.defaults),

        'call_italk' = call_italk.params %<>%
            reassign_params(., call_italk.defaults)
    )

    return(default_args)
}


#' Helper function to replace default parameters
#' @param replacements named list corresponding to default arguments
#' @param defaults named list with default arguments
#'
#' @returns a named list with the same arguments
#'
#' @noRd
reassign_params <- function(replacements,
                            defaults){
    if(!is.null(replacements)){
            nm1 <- intersect(names(defaults), names(is.na(replacements)))
            modifyList(defaults, replacements[nm1])
        } else{
            defaults
        }
}
