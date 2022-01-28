#' Function to handle different types of object as input and do basic quality checks
#'
#' @param sce Seurat or SingleCellExperiment object
#'
#' @export
liana_prep <- function (sce, ...) {
    UseMethod("liana_prep", sce)
}

#' @export
liana_prep.SingleCellExperiment <- function(sce, idents = NULL, ...){

    if(!all(c("counts", "logcounts") %in% SummarizedExperiment::assayNames(sce))){
        stop("liana expects `counts` and `logcounts` to be present in the SCE object")
    }

    idents %<>% `%||%` (SingleCellExperiment::colLabels(sce))
    if(is.null(idents)){
        stop("Please set the cell types of interest to `colLabels`")
    } else if(is.null(levels(idents))){
        idents %<>% as.factor()
        message(str_glue("`colLabels` was converted to factor"))
    }

    # Assign idents to default if not passed
    SingleCellExperiment::colLabels(sce) <- idents

    return(.filter_sce(sce))
}

#' @export
liana_prep.Seurat <- function(sce, idents = NULL, assay = NULL, ...){

    assay %<>% `%||%`(SeuratObject::DefaultAssay(sce))
    message(stringr::str_glue("Running LIANA with {assay} as default assay"))

    # Assign idents to default if not passed
    idents %<>% `%||%`(SeuratObject::Idents(sce))
    if(is.null(idents)){
        stop("Please set the cell types of interest to `Idents`")
    } else if(is.null(levels(idents))){
        idents %<>% as.factor()
        message(str_glue("`Idents` were converted to factor"))
    }

    # convert from seurat_object to sce
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(
            counts = GetAssayData(object = sce, assay = assay, slot = "counts"),
            logcounts = GetAssayData(object = sce, assay = assay, slot = "data")
        ),
        metadata = sce@meta.data)

    SingleCellExperiment::colLabels(sce) <- idents

    return(.filter_sce(sce))
}

#' Helper function to perform basic filterin on the SCE object prior to feeding it to LIANA
#'
#' @param sce SingleCellExperiment Object
#'
#' @return SingleCellExperiment object
#'
.filter_sce <- function(sce){
    # EXTEND QUALITY CONTROL STEPS
    nonzero_cells <- colSums(counts(sce)) > 0
    nonzero_genes <- rowSums(counts(sce)) > 0

    if(!all(nonzero_cells) | !all(nonzero_genes)){
        nzero_genes <- sum(map_dbl(nonzero_genes, function(x) rlang::is_false(x = x)))
        nzero_cells <- sum(map_dbl(nonzero_cells, function(x) rlang::is_false(x = x)))

        warning(
            stringr::str_glue("{nzero_genes} genes and/or {nzero_cells} ",
                              "cells were removed as they had no counts!")
        )
    }

    return(sce[nonzero_genes, nonzero_cells])
}


#' Helper function to convert sce to seurat for EXTERNAL `call_` functions only
#'
#' @param sce SingleCellExperiment or Seurat Object
#' @param assay name of the active assay
#'
#' @noRd
.liana_convert <- function(sce, assay){
    seurat_object <- SeuratObject::as.Seurat(sce)
    Idents(seurat_object) <- SingleCellExperiment::colLabels(sce)
    seurat_object@assays[[assay]] <- seurat_object@assays[[1]]
    SeuratObject::DefaultAssay(seurat_object) <- assay

    return(seurat_object)
}

