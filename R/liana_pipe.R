#' Liana Pipe which runs DE analysis and merges needed information for LR inference
#'
#' @param seurat_object
#' @param op_resource
#' @param test.type `test.type` passed to \link(scran::findMarkers)
#' @param pval.type `pval.type` passed to \link(scran::findMarkers)
#' @param seed Set Random Seed
#'
#' @import scuttle scran SingleCellExperiment Seurat
#'
#' @export
#'
#' @return Returns a tibble with information required for LR calc
liana_pipe <- function(seurat_object, # or sce object
                       op_resource,
                       test.type = "t",
                       pval.type = "all",
                       seed=1234){
    set.seed(seed)

    # Resource Format
    transmitters <- op_resource$source_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    receivers <- op_resource$target_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    entity_genes <- union(transmitters$gene, receivers$gene)


    # Filter to LigRec and scale
    seurat_object <- seurat_object[rownames(seurat_object) %in% entity_genes]
    seurat_object <- Seurat::ScaleData(seurat_object, features = entity_genes)

    # convert to SCE
    sce <- Seurat::as.SingleCellExperiment(seurat_object,  assay="RNA")
    colLabels(sce) <- Seurat::Idents(seurat_object)
    sce@assays@data$scaledata <- seurat_object@assays$RNA@scale.data

    # Get Avg Per Cluster (data assay)
    means <- scuttle::summarizeAssayByGroup(sce,
                                            ids = colLabels(sce),
                                            assay.type = "counts")
    means <- means@assays@data$mean

    # scaled (z-transformed) means
    scaled <- scuttle::summarizeAssayByGroup(sce,
                                             ids = colLabels(sce),
                                             assay.type = "scaledata")
    scaled <- scaled@assays@data$mean


    # Find Markers and Format
    cluster_markers <- scran::findMarkers(sce,
                                          groups = colLabels(sce),
                                          direction = "any",
                                          full.stats = TRUE,
                                          test.type = test.type,
                                          pval.type = pval.type,
                                          assay.type = "counts") %>%
        pluck("listData") %>%
        map(function(cluster)
            cluster %>%
                as.data.frame() %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                select(gene, p.value, FDR, stat = summary.stats))

    # Get all Possible Cluster pair combinations
    pairs <- expand_grid(source = unique(colLabels(sce)),
                         target = unique(colLabels(sce)))

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
                   entity = "ligand",
                   type = "count") %>%
        join_means(means = means,
                   source_target = "target",
                   entity = "receptor",
                   type = "count") %>%
        join_means(means = scaled,
                   source_target = "source",
                   entity = "ligand",
                   type = "scaled") %>%
        join_means(means = scaled,
                   source_target = "target",
                   entity = "receptor",
                   type = "scaled") %>%
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
#'
#' @noRd
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
#' @param type type of mean to join (count or scaled)
#'
#' @import magrittr
#'
#' @return Returns the Average Expression Per Cluster
#'
#' @noRd
join_means <- function(lr_res,
                       means,
                       source_target,
                       entity,
                       type){

    entity.avg <- sym(str_glue("{entity}.{type}"))

    means %<>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        pivot_longer(-gene, names_to = "cell", values_to = "avg") %>%
        dplyr::rename({{ source_target }} := cell,
                      {{ entity }} := gene,
                      {{ entity.avg }} := avg)

    lr_res %>%
        left_join(means, by=c(source_target, entity))
}


#' Join Expression per Cluster
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param means Gene avg expression per cluster
#' @param entity ligand or receptor
#'
#' @return Returns the Summed Average Expression Per Cluster
#'
#' @noRd
join_sum_means <- function(lr_res, means, entity){

    entity.expr = sym(str_glue("{entity}.sum"))

    sums <- means %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        rowwise("gene") %>%
        mutate(sum.means = sum(c_across(where(is.numeric)))) %>%
        select(gene, sum.means) %>%
        dplyr::rename({{ entity }} := gene,
                      {{ entity.expr }} := sum.means)  %>%
        distinct() %>%
        ungroup()

    lr_res %>%
        left_join(sums, by=c(entity))
}


