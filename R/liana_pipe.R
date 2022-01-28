#' Liana Pipe which runs DE analysis and merges needed information for LR inference
#'
#' @param sce SingleCellExperiment Object
#' @param op_resource resource tibble obtained via \link{liana::select_resource}
#' @inheritParams liana_scores
#' @inheritParams scran::findMarkers
#' @param expr_prop minimum proportion of gene expression per cell type (0.2 by default),
#'  yet one should consider setting this to an appropriate value between 0 and 1,
#'  as an assumptions of these method is that communication is coordinated at the cluster level.
#' @param assay assay to be used ("RNA" by default)
#' @param assay.type - the type of data to be used to calculate the means
#'  (counts by default), available options are: "counts" and "logcounts"
#'
#' @import SingleCellExperiment SeuratObject
#' @importFrom scran findMarkers
#' @importFrom scuttle summarizeAssayByGroup
#'
#' @export
#'
#' @return Returns a tibble with information required for LR calculations downstream
liana_pipe <- function(sce,
                       op_resource,
                       decomplexify = TRUE,
                       test.type = "wilcox",
                       pval.type = "all",
                       trim = 0,
                       assay = "RNA",
                       assay.type = "logcounts"){

    ### this whole chunk needs to move to liana_wrap
    # Resource Format
    transmitters <- op_resource$source_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    receivers <- op_resource$target_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    entity_genes = union(transmitters$gene,
                         receivers$gene)

    # calculate global_mean required for SCA
    global_mean <- Matrix::mean(
        exec(assay.type, sce)
        )

    # Filter `sce` to only include ligand receptor genes
    # and any cells which don't contain any expressed LR genes
    sce <- sce[rownames(sce) %in% entity_genes,
               Matrix::colSums(counts(sce)) > 0]
    # Scale genes across cells
    sce@assays@data[["scaledata"]] <- as.matrix(row_scale(exec(assay.type, sce)))

    # Get Avg and  Prop. Expr Per Cluster
    mean_prop <-
        scuttle::summarizeAssayByGroup(sce,
                                       ids = colLabels(sce),
                                       assay.type = assay.type,
                                       statistics = c("mean", "prop.detected"))
    means <- mean_prop@assays@data$mean
    props <- mean_prop@assays@data$prop.detected

    # scaled (z-transformed) means
    scaled <- scuttle::summarizeAssayByGroup(sce,
                                             ids = colLabels(sce),
                                             assay.type = "scaledata",
                                             statistics = c("mean"))
    scaled <- scaled@assays@data$mean

    # calculate PEM scores
    pem_scores <- compute_pem_scores(sce = sce,
                                     assay.type = assay.type)

    # Get Log2FC
    logfc_df <- get_log2FC(sce, "counts")

    # Find Markers and Format
    cluster_markers <- scran::findMarkers(sce,
                                          groups = colLabels(sce),
                                          direction = "any",
                                          full.stats = TRUE,
                                          test.type = test.type,
                                          pval.type = pval.type,
                                          assay.type = assay.type) %>%
        pluck("listData") %>%
        map2(., names(.), function(cluster, cluster_name){
            cluster %>%
                as.data.frame() %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                select(gene, p.value, FDR, stat = summary.stats)
            })

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
    lr_res %<>%
        join_means(means = means,
                   source_target = "source",
                   entity = "ligand",
                   type = "expr") %>%
        join_means(means = means,
                   source_target = "target",
                   entity = "receptor",
                   type = "expr") %>%
        join_means(means = scaled,
                   source_target = "source",
                   entity = "ligand",
                   type = "scaled") %>%
        join_means(means = props,
                   source_target = "target",
                   entity = "receptor",
                   type = "prop") %>%
        join_means(means = props,
                   source_target = "source",
                   entity = "ligand",
                   type = "prop") %>%
        join_means(means = scaled,
                   source_target = "target",
                   entity = "receptor",
                   type = "scaled") %>%
        join_sum_means(means = means,
                       entity = "ligand") %>%
        join_sum_means(means = means,
                       entity = "receptor") %>%
        # Join PEM scores
        join_means(means = pem_scores,
                   source_target = "target",
                   entity = "receptor",
                   type = "pem") %>%
        join_means(means = pem_scores,
                   source_target = "source",
                   entity = "ligand",
                   type = "pem") %>%
        # logFC
        join_log2FC(logfc_df,
                    source_target = "source",
                    entity="ligand") %>%
        join_log2FC(logfc_df,
                    source_target = "target",
                    entity="receptor") %>%
        # Global Mean
        mutate(global_mean = global_mean)

    message("LIANA: LR summary stats calculated!")

    if(decomplexify){
        # Join complexes (recomplexify) to lr_res
        cmplx <- op_resource %>%
            select(
                ligand = source_genesymbol,
                ligand.complex = source_genesymbol_complex,
                receptor = target_genesymbol,
                receptor.complex = target_genesymbol_complex
                )

        lr_res %<>%
            left_join(., cmplx,
                      by=c("ligand", "receptor")) %>%
            distinct()
    }

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
#' @importFrom magrittr %>% %<>%
#'
#' @return Returns the Average Expression Per Cluster
#'
#' @noRd
join_means <- function(lr_res,
                       means,
                       source_target,
                       entity,
                       type,
                       pb = NULL){

    if(!is.null(pb)){
        pb$tick()$print()
    }

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


#' Helper Function to join log2FC dataframe to LR_res
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param logfc_df obtained via \link{get_log2FC}
#' @param source_target target or source cell
#' @param entity ligand or receptor
#'
#' @noRd
join_log2FC <- function(lr_res,
                        logfc_df,
                        source_target,
                        entity){

    entity.fc = sym(str_glue("{entity}.log2FC"))

    logfc <- logfc_df %>%
        dplyr::rename(
            {{ source_target }} := cell,
            {{ entity }} := gene,
            {{ entity.fc }} := avg_log2FC)  %>%
        distinct()

    lr_res %>%
        left_join(logfc, by=c(entity, source_target))

}




#' Get Log2FC of Subject vs LOSO (i.e. 1 cell type vs all other cells FC)
#'
#' @param sce SingleCellExperiment object
#' @param subject leave-one-out subject, i.e. the cluster whose log2FC we wish
#'    to calculate when compared to all other cells
#' @inheritParams liana_pipe
#'
#' @return A log2FC dataframe for a given cell identity
#'
#' @details log2FC is calculated using the raw count average + a pseudocount of 1.
#' `assay.type` should be the raw counts
#'
#' @noRd
get_log2FC <- function(sce,
                       assay.type){

    # normalize counts across libraries
    sce <- scater::logNormCounts(sce,
                                 log=FALSE,
                                 assay.type=assay.type)

    # iterate over each possible cluster leaving one out
    levels(colLabels(sce)) %>%
        map(function(subject){
            # Subject (i.e. target) Cluster avg
            subject_avg <-
                scater::calculateAverage(subset(sce,
                                                select = colLabels(sce)==subject),
                                         assay.type = "normcounts"
                                         ) %>%
                as_tibble(rownames = "gene") %>%
                dplyr::rename(subject_avg = value)

            # All other cells average
            loso_avg <-
                scater::calculateAverage(subset(sce,
                                                select = !(colLabels(sce) %in% subject)),
                                         assay.type = "normcounts"
                                         ) %>%
                as_tibble(rownames = "gene") %>%
                dplyr::rename(loso_avg = value)

            # Join avg and calculate FC
            left_join(subject_avg, loso_avg, by="gene") %>%
                mutate(avg_log2FC =
                           log2((subject_avg + 1)) - log2((loso_avg + 1))) %>%
                select(gene, avg_log2FC)

        }) %>% setNames(levels(colLabels(sce))) %>%
        enframe(name = "cell") %>%
        unnest(value)
}


#' Helper Function to 'decomplexify' ligands and receptors into
#'
#' @param resource a ligrec resource
#' @param columns columns to separate and pivot long (e.g. genesymbol or uniprot)
#'
#' @return returns a longer tibble with complex subunits on seperate rows
#'
#' @noRd
decomplexify <- function(resource,
                         columns = c("source_genesymbol",
                                     "target_genesymbol")){
    columns %>%
        map(function(col){
            sep_cols <- c(str_glue("col{rep(1:5)}"))
            col.complex <- str_glue("{col}_complex")

            resource <<- resource %>%
                mutate({{ col.complex }} :=
                           resource[[str_glue("{col}")]]) %>%
                separate(col,
                         into = sep_cols,
                         sep = "_",
                         extra = "drop",
                         fill = "right") %>%
                pivot_longer(cols = all_of(sep_cols),
                             values_to = col,
                             names_to = NULL) %>%
                tidyr::drop_na(col) %>%
                distinct() %>%
                mutate_at(.vars = c(col),
                          ~str_replace(., "COMPLEX:", ""))
        })
    return(resource)
}


#' Helper Function to Get a rowwise scaled matrix
#'
#' @param mat a matrix, typically the logcounts matrix from an SCE object
#'
#' @noRd
row_scale <- function(mat){
    col_means = rowMeans(mat, na.rm = TRUE) # Get the column means
    col_sd = MatrixGenerics::rowSds(mat, center = col_means) # Get the column sd

    # return scaled mat
    return((mat - col_means) / col_sd)
}


