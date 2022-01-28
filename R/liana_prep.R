#' Function to handle different types of object
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

    # EXTEND QUALITY CONTROL STEPS

    return(sce[, colSums(counts(sce)) > 0])
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

