require(SingleCellExperiment)
require(scuttle)
require(scran)

testdata <- readRDS("inst/testdata/input/testdata.rds")

op_resource <- select_resource("OmniPath")[[1]]
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)


# Convert to SCE
test_sce <- Seurat::as.SingleCellExperiment(testdata)
colLabels(test_sce) <- Seurat::Idents(testdata)
test_sce


# Find Markers and Format
cluster_markers <- scran::findMarkers(test_sce,
                          groups = colLabels(test_sce),
                          direction = "any",
                          full.stats = TRUE,
                          test.type = "t") %>%
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


# Find Avg Per Cluster
cluster_summ <- scuttle::summarizeAssayByGroup(test_sce,
                                               ids = colLabels(test_sce))
cluster_summ@assays@data$mean

res <- pairs %>%
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

res %>%
    filter(ligand.FDR <= 0.05 & receptor.FDR <= 0.05) %>%
    mutate(stat_weight = ligand.stat * receptor.stat) %>%
    arrange(desc(stat_weight))



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
