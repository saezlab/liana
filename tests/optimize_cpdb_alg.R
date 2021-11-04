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

# Test with DT
start <- Sys.time()
dt_perms <- get_permutations_dt(pipe_out,
                                sce)
end <- Sys.time() - start
end


sce_mat <- t(as.matrix(sce@assays@data[["logcounts"]]))
col_labels <- colLabels(sce)

start <- Sys.time()
mean_agg <- stats::aggregate(sce_mat,
                             list(col_labels),
                             FUN=mean) %>%
    tibble::as_tibble() %>%
    dplyr::rename(celltype = Group.1) %>%
    tidyr::pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from=celltype,
                       id_cols=gene,
                       values_from=value) #%>%
    # tibble::column_to_rownames("gene")
end <- Sys.time() - start
end


start <- Sys.time()
xd <- scuttle::summarizeAssayByGroup(sce,
                               ids=col_labels,
                               statistics = c("mean"))
end <- Sys.time() - start
end

require(data.table)
start <- Sys.time()
dt <- data.table(value = as.vector(sce_mat),
                 gene = rownames(sce),
                 celltype = col_labels)

mean_dt <- dt[, mean(value), by=c("gene", "celltype")] %>%
    tidyr::pivot_wider(names_from=celltype,
                       id_cols=gene,
                       values_from=V1) #%>%
    # tibble::column_to_rownames("gene")
end <- Sys.time() - start
end



# convert to data table
dt <- data.table(celltype = col_labels,
                 mat = as.vector(sce_mat),
                 gene = colnames(sce_mat))
mean_dt <- dt[, mean(mat), by=c("celltype", "gene")] %>%
    tidyr::pivot_wider(names_from=celltype,
                       id_cols=gene,
                       values_from=V1) %>%
    tibble::column_to_rownames("gene")




get_permutations_dt <- function(lr_res,
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
    sce_dt <- data.table::data.table(value = as.vector(sce@assays@data[[assay.type]]),
                                     gene = rownames(sce),
                                     celltype = col_labels)

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
                       .f = mean_permute_dt,
                       sce_dt = sce_dt,
                       trim = trim,
                       pb = progress_bar,
                       parallelize = parallelize,
                       workers = workers)

    return(perm)
}


mean_permute_dt <- function(col_labels,
                            sce_dt,
                            trim,
                            pb){
    pb$tick()$print()

    sce_dt[, mean(value), by=c("celltype", "gene")] %>%
        tidyr::pivot_wider(names_from=celltype,
                           id_cols=gene,
                           values_from=V1) %>%
        tibble::column_to_rownames("gene")
}

