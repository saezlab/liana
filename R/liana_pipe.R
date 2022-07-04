#' Liana Pipe which runs DE analysis and merges needed information for LR inference
#'
#' @param sce SingleCellExperiment Object
#' @param op_resource resource tibble obtained via \link{liana::select_resource}
#' @inheritParams liana_scores
#' @inheritParams scran::findMarkers
#' @param assay assay to be used ("RNA" by default)
#' @param assay.type - the type of data to be used to calculate the means
#'  (logcounts by default), available options are: "counts" and "logcounts"
#' @param verbose logical for verbosity
#' @param cell.adj cell adjacency tibble/dataframe /w weights by which we will
#' `multiply` the relevant columns. Any cell pairs with a weights of 0 will be
#' filtered out.
#' Note that if working with LIANA's default methods, we suggest weights >= 0 & =< 1.
#' This ensure that all methods' score will be meaningfully weighed without
#' changing the interpretation of their scores, thus allow one to filter SCA,
#' rank NATMI, etc.
#'
#' @inheritParams .antilog1m
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
                       test.type = "wilcox",
                       pval.type = "all",
                       assay = "RNA",
                       assay.type = "logcounts",
                       verbose = TRUE,
                       base,
                       cell.adj = NULL){

    ### this whole chunk needs to move to liana_wrap
    # Resource Format
    transmitters <- op_resource$source_genesymbol %>%
        as_tibble() %>%
        select(gene = value)
    receivers <- op_resource$target_genesymbol %>%
        as_tibble() %>%
        select(gene = value)

    # calculate global_mean required for SCA
    global_mean <- fast_mean(exec(assay.type, sce))

    # Filter `sce` to only include ligand receptor genes
    # and exclude any cells with 0 counts of LR genes
    sce <- .prep_universe(sce,
                          entity_genes = union(transmitters$gene,
                                               receivers$gene),
                          verbose)

    # Get Log2FC (done after non-expr. cell and gene filter from `liana_prep`)
    # also any cells with 0 counts of LR genes are removed (in`.prep_universe`)
    logfc_df <- get_log2FC(
        sce,
        assay.type = assay.type,
        base
    )


    # Scale genes across cells
    sce@assays@data[["scaledata"]] <- row_scale(exec(assay.type, sce))

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

    liana_message("LIANA: LR summary stats calculated!",
                  verbose = verbose
    )

    # Weigh by (spatial) constrains
    if(!is.null(cell.adj)){
        lr_res %<>%
            .sp_costrain(cell.adj)
    }


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
#' @param pb progress bar
#'
#' @importFrom magrittr %>% %<>%
#'
#' @return Returns the Average Expression Per Cluster
#'
#' @export
#'
#' @keywords internal
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

    entity.expr <- sym(str_glue("{entity}.sum"))

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

    entity.fc <- sym(str_glue("{entity}.log2FC"))

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
#' @inheritParams .antilog1m
#'
#' @inheritParams liana_pipe
#'
#' @return A log2FC dataframe for a given cell identity
#'
#' @details log2FC is calculated using the raw count average + a pseudocount of 1.
#' `assay.type` should be the raw counts
#'
#' @noRd
get_log2FC <- function(sce, assay.type, base){


    if(!is.nan(base)){
        # Get anti-logged normalized counts (preserves batch effets)
        sce@assays@data[["normcounts"]] <- .get_invcounts(sce, assay.type, base)
    } else{
        # Get raw counts as they are
        sce@assays@data[["normcounts"]] <- counts(sce)
    }


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


#' Helper Function to 'decomplexify' ligands and receptors into individual subunits
#'
#' @param resource a ligrec resource
#'
#' @param columns columns to separate and pivot long (e.g. genesymbol or uniprot),
#' `source_genesymbol` and `target_genesymbol` by default
#'
#' @return returns a longer tibble with complex subunits on seperate rows
#'
#' @details takes any number of columns, and assumes `_` as sep.
#'
#' @export
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
                tidyr::drop_na(all_of(col)) %>%
                distinct() %>%
                mutate(across(all_of(c(col, col.complex)),
                              ~str_replace(., "COMPLEX:", "")))
        })
    return(resource)
}


#' Helper Function to Get a rowwise scaled matrix
#'
#' @param mat a matrix, typically the logcounts matrix from an SCE object
#'
#'
#' @noRd
row_scale <- function(mat){
    col_means = rowMeans(mat,
                         na.rm = TRUE) # Get the column means
    col_sd = MatrixGenerics::rowSds(mat,
                                    center = col_means,
                                    na.rm = TRUE) # Get the column sd

    # return scaled mat
    return(as.matrix((mat - col_means) / col_sd))
}

#' Helper function to inverse logged counts
#'
#' @param x mat or array
#' @param base base for conversion from log-tranformed ~CPM back to ~CPM.
#'
#' @keywords internal
.antilog1m <- function(x, base=2){base ^ (x) - 1}


#' Helper function to generate inversed counts (i.e. normalized but not logged)
#'
#' @param sce SingleCellExperiment object
#' @param assay.type counts slot
#' @param base a positive or complex number: the base with respect to which
#' log-transformation was computed.
#'
#' @noRd
.get_invcounts <- function(sce, assay.type, base){
    antilogged <- .antilog1m(slot(exec(assay.type, sce), "x"), base = base)

    methods::new(
        "dgCMatrix",
        i = slot(exec(assay.type, sce), "i"),
        p = slot(exec(assay.type, sce), "p"),
        Dim = slot(exec(assay.type, sce), "Dim"),
        Dimnames = slot(exec(assay.type, sce), "Dimnames"),
        x = antilogged,
        factors = slot(exec(assay.type, sce), "factors")
    )
}


#' Helper function to fix r mean
#'
#' @param mat a matrix
#'
#' @details r mean is slow and it overflows on memory.
#'
#' @noRd
fast_mean <- function(mat){
    if(class(mat)=="dgCMatrix"){
        sum(mat@x)/(as.numeric(nrow(mat)) * as.numeric(ncol(mat)))
    } else{
        Matrix::mean(mat)
    }
}

#' Helper function to format SCE to the LR universe
#' @param sce SingleCellExperiment object
#' @param entity_genes union of ligand-receptor genes
#' @param verbose verbose - True/Flase
#'
#' @noRd
.prep_universe <- function(sce, entity_genes, verbose){

    # Keep only LR universe
    sce <- sce[rownames(sce) %in% entity_genes, ]

    # Check organism/gene intersect
    if(nrow(sce) < 3){
        liana_message(
            "Low gene intersect (<3) detected!",
            "Please check if the rownames of the data match the gene identities in the resource (i.e. human genesymbols).",
            output = "stop",
            verbose = verbose
        )
    }

    # Check for non-zero cells
    nonzero_cells <- colSums(counts(sce)) > 0

    if(!all(nonzero_cells)){
        nzero_cells <- sum(map_dbl(nonzero_cells, function(x) rlang::is_false(x = x)))

        liana_message(
            stringr::str_glue("{nzero_cells} cells were excluded as they",
                              " did not express any ligand-receptor genes!"),
            output="warning",
            verbose=verbose
        )
    }

    return(sce[,nonzero_cells])
}


#' Weigh by spatial constrans
#'
#' @param lr_res liana_pipe output prior to joining complexes
#' @param cell.adj cell adjacency weights (should be positive)
#' @param adjacency name of the column with cell pair adjacency scores
#'
#' @return weighed lr_res
#'
#' @details Note that for the case that there are weights from 0-1, the negative
#' values (in e.g. logFC, z-scores) might be counter-logically affected - i.e.
#' they would be brought closer to 0.
#' Thus, by default liana expects weights from 0-1. These are then multiplied
#' for positive values, while negative values are divided.
#'
#' Alternatively, one could e.g. multiply the weights by a factor (e.g. 10,000),
#' if logFC and Connectome are used. However, this would change some of
#' the assumptions/interpretations of the scores, while
#' consensus ranking will be unaffected.
#'
#' Also, note that any interactions between cell pairs with an adjacency of 0
#' will be excluded (this would affect the scores from CytoTalk).
#'
#' NB! `%/*/%` is only applicable and relevant to logFC and z-scores from a
#' single-context, and should not be used when scaling between conditions!!!
#'
#' @keywords internal
.sp_costrain <- function(lr_res,
                         cell.adj,
                         adjacency = "adjacency"){
    lr_res %>%
        left_join(cell.adj, by=c("source", "target")) %>%
        # remove any interactions between non-interacting cells
        filter(.data[[adjacency]]!=0) %>%
        # weigh relevant columns by adjacency
        mutate(across(ends_with(c("expr", "pem") #**
        ), ~`*`(.x, adjacency))) %>%
        rowwise() %>% # required as its done by element (not vector operation)
        mutate(across(ends_with(c("log2FC", "scaled")),
                      # makes sense for single context
                      # but makes no sense for multiple (better to change to prod)
                      ~`%/*/%`(.x, adjacency))) %>%
        ungroup()
}



#' Helper function to multiply or divide depending on sign
#'
#' @param x target for weighing
#' @param y weight
#'
#' @noRd
`%/*/%` <- function(x, y){
    if(y==0) return(0)

    if(x >= 0){
        x * y
    } else{
        x / y
    }
}
