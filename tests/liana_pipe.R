#' Liana Pipe function that returns information required for LR calc
#'
#' @param seurat_object
#' @param op_resource
#' @param test.type `test.type` passed to \link(scran::findMarkers)
#'
liana_pipe <- function(seurat_object, # or sce object
                       op_resource,
                       test.type = "t"){

    # Convert to SCE
    test_sce <- Seurat::as.SingleCellExperiment(seurat_object,
                                                assay="RNA")
    colLabels(test_sce) <- Seurat::Idents(seurat_object)
    # test_sce <- scuttle::logNormCounts(test_sce)

    # Format OmniPath
    transmitters <- op_resource$source_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    receivers <- op_resource$target_genesymbol %>%
        as_tibble() %>%
        select(gene = value)

    # Find Markers and Format
    cluster_markers <- scran::findMarkers(test_sce,
                                          groups = colLabels(test_sce),
                                          direction = "any",
                                          full.stats = TRUE,
                                          test.type = test.type) %>%
        pluck("listData") %>%
        map(function(cluster)
            cluster %>%
                as.data.frame() %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                select(gene, p.value, FDR, stat = summary.stats))


    # Get all Possible Cluster pair combinations
    pairs <- expand_grid(source = unique(colLabels(test_sce)),
                         target = unique(colLabels(test_sce)))

    # Get Avg Per Cluster (data assay)
    means <- scuttle::summarizeAssayByGroup(test_sce,
                                            ids = colLabels(test_sce),
                                            assay.type = "counts")
    means <- means@assays@data$mean


    # Get DEGs to LR format
    lr_res <- pairs %>%
        pmap(function(source, target){
            source_stats <- ligrec_degformat(cluster_markers[[source]],
                                             entity = transmitters,
                                             source_target = "source")
            target_stats <- ligrec_degformat(cluster_markers[[target]],
                                             entity = receivers,
                                             source_target = "target")

            op_resource %>%
                select(ligand = source_genesymbol,
                       receptor = target_genesymbol) %>%
                left_join(source_stats, by = "ligand") %>%
                left_join(target_stats, by = "receptor") %>%
                na.omit() %>%
                distinct() %>%
                mutate(source = source,
                       target = target)
        }) %>%
        bind_rows()

    # Join Expression Means
    lr_res <- lr_res %>%
        join_means(means = means,
                   source_target = "source",
                   entity = "ligand") %>%
        join_means(means = means,
                   source_target = "target",
                   entity = "receptor") %>%
        select(source, starts_with("ligand"),
               target, starts_with("receptor"),
               everything()) %>%
        join_sum_means(means = means,
                       entity = "ligand") %>%
        join_sum_means(means = means,
                       entity = "receptor") %>%
        select(source, starts_with("ligand"),
               target, starts_with("receptor"),
               everything())

    return(lr_res)
}



#' Helper Function to join DEG stats to LR
#'
#' @param cluster_markers dataframe with DE stats for a cluster
#' @param entity Transmitter or Receiver vector passed as tibble
#' @param source_target whether this is the source or target cluster
#'
#' @return A tibble with stats for receivers or transmitters per cluster
ligrec_degformat <- function(cluster_markers,
                             entity,
                             source_target){
    cluster_markers %>%
        left_join(entity, ., by = "gene") %>%
        na.omit() %>%
        {
            if(source_target=="source"){
                dplyr::select(
                    .,
                    ligand = gene,
                    ligand.pval = p.value,
                    ligand.FDR = FDR,
                    ligand.stat = stat
                )
            }else if(source_target=="target"){
                dplyr::select(
                    .,
                    receptor = gene,
                    receptor.pval = p.value,
                    receptor.FDR = FDR,
                    receptor.stat = stat
                )
            } else{
                stop("Incorrect entity!")
            }
        } %>%
        distinct() %>%
        na.omit()
}




#' Join Expression per Cluster
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param means Gene avg expression per cluster
#' @param source_target target or source cell
#' @param entity ligand or receptor
#'
#' @return Returns the Average Expression Per Cluster
join_means <- function(lr_res, means, source_target, entity){

    entity.avg = sym(str_glue("{entity}.avg"))

    means_pivot <- means %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        pivot_longer(-gene, names_to = "cell", values_to = "avg") %>%
        dplyr::rename({{ source_target }} := cell,
                      {{ entity }} := gene,
                      {{ entity.avg }} := avg)

    lr_res %>%
        left_join(means_pivot, by=c(source_target, entity))
}


#' Join Expression per Cluster
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param means Gene avg expression per cluster
#' @param entity ligand or receptor
#'
#' @return Returns the Summed Average Expression Per Cluster
join_sum_means <- function(lr_res, means, entity){

    # Sum of the mean expression across all cell types
    entity.expr = sym(str_glue("{entity}.sum"))

    sums <- means %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        rowwise("gene") %>%
        mutate(sum.means = sum(c_across(where(is.numeric)))) %>%
        select(gene, sum.means) %>%
        dplyr::rename({{ entity }} := gene,
                      {{ entity.expr }} := sum.means)  %>%
        distinct()

    lr_res %>%
        left_join(sums, by=c(entity))
}


