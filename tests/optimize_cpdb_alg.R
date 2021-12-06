# load liana pipe out
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
pipe_out <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_pipe.RDS"))
op_resource <- select_resource("OmniPath")[[1]]

sce <- seurat_to_sce(seurat_object,
                     entity_genes = union(op_resource$target_genesymbol,
                                          op_resource$source_genesymbol),
                     assay="RNA")

# Test with aggregate
start <- Sys.time()
agg_perms <- get_permutations(pipe_out,
                              sce)
end <- Sys.time() - start
end

# Test with SCE
start <- Sys.time()
sce_perms <- get_permutations_sce(pipe_out,
                                  sce)
end <- Sys.time() - start
end


sce_mat <- t(as.matrix(sce@assays@data[["logcounts"]]))

start <- Sys.time()
mean_agg <- stats::aggregate(sce_mat,
                             list(colLabels(sce)),
                             FUN=mean) %>%
    tibble::as_tibble() %>%
    dplyr::rename(celltype = Group.1) %>%
    tidyr::pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from=celltype,
                       id_cols=gene,
                       values_from=value) %>%
    arrange("gene", "B", "CD8 T", "NK") # %>%
    # tibble::column_to_rownames("gene")
end <- Sys.time() - start
end


start <- Sys.time()
mean_sce <- scuttle::summarizeAssayByGroup(sce,
                                           ids=colLabels(sce),
                                           statistics = c("mean"),
                                           assay.type = "logcounts")@assays@data$mean %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble() %>%
    arrange("gene", "B", "CD8 T", "NK")
end <- Sys.time() - start
end



# sce

#' Helper Function to generate shuffled means
#'
#' @param lr_res liana_pipe results
#' @param sce SingleCellExperiment Object
#' @param nperms number of permutations
#' @param seed number used to set random seed
#' @inheritParams liana_pipe
#' @inheritParams map_custom
#'
#' @return Returns a list of shuffled gene means by cluster
#'
#' @details This function could be made generalizable to any set of genes,
#'   depending on the set (currently lr_res genes) that is used to filter - i.e.
#'   it could be replaced with e.g. genes from TF regulons
get_permutations <- function(lr_res,
                             sce,
                             nperms = 1000,
                             seed = 1234,
                             trim = 0.1,
                             parallelize = FALSE,
                             workers = 4,
                             assay.type = "logcounts"){
    # remove genes absent in lr_res
    lr_genes <- union(lr_res$ligand, lr_res$receptor)
    sce <- sce[rownames(sce) %in% lr_genes, ]

    # shuffle columns
    set.seed(seed)
    shuffled_clusts <- map(1:nperms, function(perm){
        colLabels(sce) %>%
            as_tibble(rownames = "cell") %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })

    # progress_bar
    progress_bar <- progress_estimated(nperms)

    # generate mean permutations
    perm <- map_custom(.x = shuffled_clusts,
                       .f = mean_permute_sce,
                       sce = sce,
                       trim = trim,
                       pb = progress_bar,
                       parallelize = parallelize,
                       workers = workers)

    return(perm)
}


#' Function to calculate mean LR expression from shuffled cluster label matrices
#'  as done in CellPhoneDB
#'
#' @param sce_matrix single cell expression matrix (transposed)
#' @param col_labels cluster labels
#' @param pb progress bar object
#'
#' @importFrom dplyr progress_estimated
#'
#' @return Returns a list of means per gene calculated with reshuffled
#'    cluster/cell identity labels
mean_permute <- function(col_labels,
                         sce,
                         pb){
    pb$tick()$print()

    scuttle::summarizeAssayByGroup(sce,
                                   ids=col_labels,
                                   statistics = c("mean"),
                                   assay.type = "logcounts")@assays@data$mean
}



