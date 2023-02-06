#' Function to handle different types of object as input and do basic quality checks
#'
#' @param sce SingleCellExperiment or Seurat object
#' @param ... dot dot dot bucket - not passed to anything, handles issues with
#' passing non-existing arguments
#'
#' @export
liana_prep <- function (sce, ...) {
    UseMethod("liana_prep", sce)
}

#' @export
liana_prep.SingleCellExperiment <- function(sce,
                                            idents_col = NULL,
                                            verbose = TRUE,
                                            min_cells = 0,
                                            ...){

    if(!all(c("counts", "logcounts") %in% SummarizedExperiment::assayNames(sce))){
        stop("liana expects `counts` and `logcounts` to be present in the SCE object")
    }

    # Assign idents to default if not passed
    idents <-
        .format_idents(metadata = as.data.frame(SingleCellExperiment::colData(sce)),
                       active_idents = SingleCellExperiment::colLabels(sce),
                       idents_col = idents_col,
                       object_class = class(sce) %>% pluck(1),
                       verbose = verbose)

    # Assign idents to default if not passed
    SingleCellExperiment::colLabels(sce) <- idents
    sce@int_metadata$base <- 2 # save base for logFC conv

    return(.filter_sce(sce, min_cells, verbose))
}

#' @export
liana_prep.Seurat <- function(sce,
                              idents_col = NULL,
                              verbose = TRUE,
                              assay = NULL,
                              min_cells = 0,
                              ...){

    assay %<>% `%||%`(SeuratObject::DefaultAssay(sce))
    message(stringr::str_glue("Expression from the `{assay}` assay will be used"))

    # Assign idents to default if not passed
    # Assign idents to default if not passed
    idents <-
        .format_idents(metadata = sce@meta.data,
                       active_idents = SeuratObject::Idents(sce),
                       idents_col = idents_col,
                       object_class = class(sce) %>% pluck(1),
                       verbose = verbose)

    # convert from seurat_object to sce
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(
            counts = SeuratObject::GetAssayData(object = sce,
                                                assay = assay,
                                                slot = "counts"),
            logcounts = SeuratObject::GetAssayData(object = sce,
                                                   assay = assay,
                                                   slot = "data")
        ),
        metadata = sce@meta.data)

    SingleCellExperiment::colLabels(sce) <- idents
    sce@int_metadata$base <- exp(1) # save base for logFC conv

    return(.filter_sce(sce, min_cells, verbose))
}


#' Helper function to perform basic filterin on the SCE object prior to feeding it to LIANA
#'
#' @param sce SingleCellExperiment Object
#' @param min_cells minimum cell per cell identity to be considered for analysis
#' @param verbose logical for verbosity
#'
#' @return SingleCellExperiment object
#'
#' @noRd
.filter_sce <- function(sce, min_cells, verbose){
    # EXTEND QUALITY CONTROL STEPS
    if(any(is.na(colLabels(sce)))){
        stop("NAs found in Idents/Labels!")
    }

    if (!is.numeric(min_cells) | min_cells <= 0) {
      stop("min_cells parameter should be a numeric value greater than 0")
    }
    
    # Remove any cell types with less than X cells
    remove_labels <- colData(sce) %>%
        as_tibble() %>%
        mutate(label = as.character(label)) %>%
        group_by(label) %>%
        summarise(remove=n() < min_cells) %>%
        filter(remove) %>%
        pull(label)
    if(length(remove_labels) > 0){
        liana_message(
            stringr::str_glue(
                "Cell identities with less ",
                "than {min_cells} cells: {remove_labels} were removed!"
                ),
            output="message",
            verbose=verbose
        )
        sce <- sce[,!colLabels(sce) %in% remove_labels]
        colLabels(sce) <- as.factor(as.character(colLabels(sce)))
    }

    nonzero_genes <- rowSums(counts(sce)) > 0
    nonzero_cells <- colSums(counts(sce)) > 0

    if(!all(nonzero_cells) | !all(nonzero_genes)){
        nzero_genes <- sum(map_dbl(nonzero_genes, function(x) rlang::is_false(x = x)))
        nzero_cells <- sum(map_dbl(nonzero_cells, function(x) rlang::is_false(x = x)))

        liana_message(
            stringr::str_glue("{nzero_genes} genes and/or {nzero_cells} ",
                              "cells were removed as they had no counts!"),
            output="warning",
            verbose=verbose
        )
    }

    # Convert to sparse!
    if(!is(counts(sce), "sparseMatrix")){
        counts(sce) <- as(counts(sce), "sparseMatrix")
        logcounts(sce) <- as(logcounts(sce), "sparseMatrix")
    }

    # Check counts (if negative -> Stop)
    if(min(exec("logcounts", sce)) < 0){
        stop("Negative counts are present in the Matrix!")
    }

    # Check counts (if not log-transformed -> Stop)
    if(all(round(exec("logcounts", sce)@x[1:100]) == exec("logcounts", sce)@x[1:100])){
        stop("Please make sure an assay with normalized counts is present!")
    }


    return(sce[nonzero_genes, nonzero_cells])
}



#' Helper Function to get/format the required indentity if required
#'
#' @param metadata metadata obtained from the sce object
#' @param active_idents active idents (e.g. Seurat::Idents, or sce colLabels)
#' @param object_class class of the sce object (SingleCellExperiment or Seurat)
#'
#' @inheritParams liana_wrap
#'
#' @noRd
.format_idents <- function (metadata,
                            active_idents,
                            idents_col,
                            object_class,
                            verbose){
    if (is_null(idents_col)){
        # get active ident col and assign
        idents_col <- .get_ident(metadata,
                                 active_idents,
                                 object_class = object_class)
    }

    if(!is_null(idents_col)){
        # If idents_col is not null set it that one
        if(!idents_col %in% colnames(metadata)){
            stop(str_glue("`{idents_col}` was not found!"))
        }
        idents <- metadata[[idents_col]]
        liana_message(str_glue("Running LIANA with `{idents_col}` as labels!"),
                      verbose = verbose)

    } else if(!is_null(active_idents)){
        idents <- active_idents
        liana_message(str_glue("Running LIANA with `colLabels`/`Idents` as labels",
                               " (matching column in metadata not found).",
                               .sep = ""),
                      verbose = verbose)
    } else{
        stop("Please provide existing cell type identities!")
    }

    # Check if idents is a factor
    if(is.null(levels(idents))){
        idents %<>% as.factor()
        message(str_glue("`Idents` were converted to factor"))

    }


    return(idents)
}




#' Helper Function to get the required identity (cluster annotation column)
#'  from the sce/seurat object metadata
#'
#' @param metadata df/tibble with metadata information
#'
#' @noRd
.get_ident <- function(metadata, idents, object_class){
    map(names(metadata),
        function(x){
            p <- metadata %>%
                select(sym(x)) %>%
                {`if`(object_class=="SingleCellExperiment",
                     .,
                     rownames_to_column(., "names")
                     )} %>%
                deframe()

            if(identical(p, idents)){
                return(x)
            }
            return()
        }) %>%
        compact %>%
        as.character %>%
        pluck(1) # to handle scenario when there are two identical columns
}



#' LIANA message/warning helper function to allow for verbosity
#'
#' @inheritParams base::stop
#' @param output type of output - message, warning, or stop
#' @param verbose logical for verbosity
liana_message <- function(...,
                          output = "message",
                          verbose = TRUE){
    if(verbose){
        exec(output, ...)
    }
}



#' Helper function to convert sce to seurat for EXTERNAL `call_` functions only
#'
#' @param sce SingleCellExperiment or Seurat Object
#' @param assay name of the active assay
#'
#' @noRd
.liana_convert <- function(sce, assay){
    seurat_object <- SeuratObject::as.Seurat(sce)
    SeuratObject::Idents(seurat_object) <- SingleCellExperiment::colLabels(sce)
    seurat_object@assays[[assay]] <- seurat_object@assays[[1]]
    SeuratObject::DefaultAssay(seurat_object) <- assay

    return(seurat_object)
}
