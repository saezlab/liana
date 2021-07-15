
# Reimplement italk FindLR
idents <- as.character(unique(Idents(seurat_object)))
pairs <- expand_grid(source = idents, target = idents)

# pairs <- combn(idents, 2) %>% t %>% as_tibble() %>% rename(source=V1, target=V2)


op_resource <- select_resource("OmniPath")[[1]]
# transmitter or receiver (iTALK reimplementation)
transmitter <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receiver <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)

res2 <- pairs %>%
    pmap(function(source, target){
        source <- format_cell(deg[[source]], transmitter, "source")
        target <- format_cell(deg[[target]], receiver, "target")

        op_resource %>%
            select(ligand = source_genesymbol,
                   receptor = target_genesymbol) %>%
            left_join(source, by = "ligand") %>%
            left_join(target, by = "receptor") %>%
            na.omit() %>%
            distinct()
    }) %>%
    bind_rows()







# Loop over all pairs
# find_lr <- function()
# Args:
# Source DE, Target DE, op_resource, ... for Seurat::FindAllMarkers
# Transmiters and Receivers info from DEG
# Target and Source LR
# full join
# Format


# Helper function to Format Target and Source
format_cell <- function(cluster_markers,
                        entity,
                        source_target){
    cluster_markers %>%
        select(
            p.value = 'p_val',
            logFC = 'avg_logFC',
            q.value = 'p_val_adj',
            cell_type = 'cluster',
            gene = 'gene'
        ) %>%
        left_join(entity, ., by = "gene") %>%
        na.omit() %>%
        {
            if(source_target=="source"){
                dplyr::select(
                    .,
                    ligand = gene,
                    cell_from_logFC = logFC,
                    cell_from_q.value = q.value,
                    cell_from = cell_type
                )
            }else if(source_target=="target"){
                dplyr::select(
                    .,
                    receptor = gene,
                    cell_to_logFC = logFC,
                    cell_to_q.value = q.value,
                    cell_to = cell_type
                )
            } else{
                stop("Incorrect entity!")
            }
        }
}



