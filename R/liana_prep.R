#' Function to handle different types of object
#' @param sce Seurat or SingleCellExperiment object
#'
#' @export
liana_prep <- function (sce, ...) {
    UseMethod("liana_prep", sce)
}

#' @export
liana_prep.SingleCellExperiment <- function(sce, identity = NULL, ...){
    # Assign identity to default if not passed
    identity %<>% `%||%`(SingleCellExperiment::colLabels(sce))
    SingleCellExperiment::colLabels(sce) <- identity

    # EXTEND QUALITY CONTROL STEPS

    return(sce[, colSums(counts(sce)) > 0])
}

#' @export
liana_prep.Seurat <- function(sce, identity = NULL, ...){
    # Assign identity to default if not passed
    identity %<>% `%||%`(Idents(sce))

    # convert from seurat_object to sce
    sce <- Seurat::as.SingleCellExperiment(sce)
    SingleCellExperiment::colLabels(sce) <- identity

    # EXTEND QUALITY CONTROL STEPS

    return(sce[, colSums(counts(sce)) > 0])
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

