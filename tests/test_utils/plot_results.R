#' Liana dotplot interactions by source and target cells
#'
#' @param liana_agg aggregated `liana_wrap` results -> preferentially filtered
#' by some condition (e.g. preferential ranking, specific interactions, etc)
#'
#' @param source_groups names of the source(sender) cell types
#' @param target_groups names of the target cell types
#'
#' @param specificity column to represent the specificity of the interaction
#' @param magnitude column to represent the magnitude of interaction (by default
#' 'sca.LRscore')
#'
#' @details Here, we refer to `specificity` as how specific this interaction is
#' to a cell type pair regards to the rest of the cell type pairs (
#' e.g. CellPhoneDB's p-values, NATMI's specificity edges, Connectome's scaled weights, etc)
#'
#' `magnitude` on the other hand is a direct measure of the expression alone,
#' by default we use SingleCellSignalR's dataset indepent LRscore (bound between 0 and 1).
#' Yet, one could also use CellChat's probabilities or CellPhoneDB's means, etc.
#'
#' @import ggplot2
liana_dotplot <- function(liana_agg,
                          source_groups,
                          target_groups,
                          specificity = "natmi.edge_specificity",
                          magnitude = "sca.LRscore"){

    # Modify for the plot
    liana_mod <- liana_agg %>%
        # Filter to only the cells of interest
        filter(source %in% source_groups) %>%
        filter(target %in% target_groups) %>%
        rename(magnitude = !!magnitude) %>%
        rename(specificity = !!specificity) %>%
        unite(c("ligand", "receptor"), col = "interaction", sep = " -> ") %>%
        unite(c("source", "target"), col = "source_target", remove = FALSE)

    # colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    cbPalette <- c("#E69F00", "#56B4E9",
                   "#009E73", "#F0E442", "#0072B2",
                   "#D55E00", "#CC79A7")

    # plot
    suppressWarnings(
        ggplot(liana_mod,
               aes(x = interaction,
                   y = target,
                   colour = specificity,
                   size = magnitude,
                   group = target
               )) +
            geom_point() +
            scale_color_gradientn(colours = viridis::viridis(20)) +
            scale_size_continuous(range = c(5, 9)) +
            facet_grid(source ~ .,
                       space = "free",
                       scales ="free",
                       switch="y") +
            theme_bw(base_size = 20) +
            theme(
                legend.text = element_text(size = 16),
                axis.text.y = element_text(colour =
                                               cbPalette[1:length(
                                                   unique(liana_mod$source)
                                               )],
                                           face = "bold",
                                           size = 23),
                axis.text.x = element_text(size = 18,
                                           angle = 90,
                                           vjust = 0.5),
                legend.title = element_text(size = 18),
                panel.spacing = unit(0.1, "lines"),
                strip.background = element_rect(fill = NA),
                strip.text = element_text(size = 24, colour = "gray6") #,
                # strip.text.y.left = element_text(angle = 0)
            ) +
            scale_y_discrete(position = "right") +
            labs(x = "Interactions (Ligand -> Receptor)",
                 colour = "Expression\nMagnitude",
                 size = "Interaction\nSpecificity",
                 y = NULL
            )
    )
}



#' Function to get filtered, processed, and log-transform pseudobulk counts
#'
#'  @param seurat_object Seurat object with celltypes
#'  @param assay Assay to consider
#'  @param expr_prop minimum proportion of gene expression per cell type
#'  @param sum_count_thresh minimum (summed) counts per gene
#'
#'  @return returns a tibble with nested counts and logcounts per celltype
get_pseudobulk <- function(seurat_object,
                           assay,
                           expr_prop,
                           sum_count_thresh){
    # show celltypes considered
    levels(Idents(seurat_object)) %>%
        map(function(lev) message(lev))

    # Convert Seurat Object to SCE
    sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)
    colLabels(sce) <- SeuratObject::Idents(seurat_object)
    gc()

    # Pseudobulk - sum all counts by celltype (+ gene expression proportions)
    pseudobulk <- scuttle::summarizeAssayByGroup(sce,
                                                 ids = colLabels(sce),
                                                 assay.type = "counts", # raw counts
                                                 statistics = c("sum", "prop"))
    pseudobulk_expr <- pseudobulk@assays@data$sum %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "sum_count")
    pseudobulk_prop <- pseudobulk@assays@data$prop.detected %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "prop")

    # Filter according to expression
    pseudo <- pseudobulk_expr %>%
        left_join(pseudobulk_prop, by = c("gene", "celltype")) %>%
        # filter genes not expressed in at least 10% of cells per celltype
        # and only keep those with summed count of at least 10
        filter(prop >= expr_prop) %>%
        filter(sum_count >= sum_count_thresh) %>%
        select(-prop)

    # looks good
    pseudo %>%
        mutate(celltype = as.factor(celltype)) %>%
        ggplot(aes(x=celltype, y = log2(sum_count))) +
        geom_violin(trim=FALSE)

    # Nest by Celltype, format, and normalize
    pseudo %<>%
        group_by(celltype) %>%
        group_nest(.key = "counts") %>%
        # format and log transform
        mutate(logcounts = counts %>%
                   map(function(c) c %>%
                           as.data.frame() %>%
                           column_to_rownames("gene") %>%
                           as.matrix() %>%
                           log2()))
}

#' Function to load and format cytosig signatures
#'
#' @param cytosig_path path to cytosig matrix
#'
#' @details reads cytosig matrix and converts it to long tibble in decoupleR
#' network format
load_cytosig <- function(cytosig_path = "data/input/cytosig/cytosig_signature_centroid.csv",
                         n_gene = 500){
    ## read CytoSig
    # models/signatures were obtained from https://github.com/data2intelligence/CytoSig/tree/master/CytoSig
    cyto_signatures <- read.table(cytosig_path,
                                  header = TRUE,
                                  row.names = 1)

    # Format Cytosig
    cytosig_net <- cyto_signatures %>%
        as.data.frame() %>%
        rownames_to_column("target") %>%
        pivot_longer(-target,
                     names_to = "cytokine",
                     values_to = "weight") %>%
        # keep top 500 genes per cytokine
        group_by(cytokine) %>%
        slice_max(n=n_gene, order_by = abs(weight)) %>%
        select(cytokine, target, weight) %>%
        mutate(mor = if_else(weight>=0, 1, -1)) %>%
        mutate(weight = abs(weight)) %>%
        mutate(cytokine = if_else(cytokine=="A",
                                  "Activin_A",
                                  cytokine)) %>%
        ungroup()

    return(cytosig_net)
}


# Cytosig: add aliases (as in OP) and cytokine family members if appropriate
aliases_and_families <- list("CD40L" = "CD40LG",
                             "GSFC" = "CSF3",
                             "IFN1" = c("IFNA1", "IFNA2", "IFNA10",
                                        "IFNA7", "IFNA21", "IFNA5",
                                        "IFNA14", "IFNA17", "IFNA6",
                                        "IFNA4",  "IFNA16", "IFNA8" ),
                             "IFNL" = c("IFNL1", "IFNL2", "IFNL3", "IFNL4"),
                             "IL12" = c("IL12A", "IL12B"),
                             "IL36" = c("IL36A", "IL36B", "IL36G", "IL36RN"),
                             "MCSF" = "CSF1",
                             "TNFA" = c("TNF", "TNFA"),
                             "TRAIL" = c("TNFSF10"),
                             "TWEAK" = "TNFSF12")
alias_tib <- tibble(cytokine = names(aliases_and_families),
                    aliases = aliases_and_families) %>%
    unnest(aliases)





# Install ggplot2 if not present
if(!requireNamespace("ggplot2")){
    install.packages("ggplot2")
}
library(ggplot2)
# devtools::install_github("saezlab/decoupleR")
require(decoupleR)

# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))

# Run liana
liana_res <- liana_wrap(seurat_object)
liana_res %<>% liana_aggregate()

# Load Cytosig
cytosig_net <- load_cytosig("~/Repos/ligrec_decouple/data/input/cytosig/cytosig_signature_centroid_expand.csv",
                            n_gene = 1000)

# Get pseudobulk at the celltype level
pseudo <- get_pseudobulk(seurat_object,
                         assay = "RNA",
                         expr_prop = 0.1,
                         sum_count_thresh = 5)

# Cytokine Activity Enrichment
pseudo_cytosig <- pseudo %>%
    mutate(cytosig_res = logcounts %>%
               map(function(logc){
                   run_mlm(
                       logc,
                       cytosig_net,
                       .source = "cytokine",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "weight",
                       sparse = FALSE
                   ) %>%
                       # rename
                       select(cytokine=source,
                              NES=score,
                              p_value)
               })
    ) %>%
    select(celltype, cytosig_res) %>%
    unnest(cytosig_res) %>%
    mutate(p_adj = p.adjust(p_value)) %>%
    arrange(p_adj)

# Obtain Active Cytokines in given Celltypes
active_cytokines <- pseudo_cytosig %>%
    filter(NES > 0 # & p_adj < 0.05 # uncomment this, to filter p_adj < 0.05 (I use test data)
           ) %>%
    print() %>%
    select(celltype, cytokine) %>%
    mutate(active_flag = 1) %>%
    left_join(alias_tib) %>%
    mutate(aliases = if_else(is.na(aliases),
                             cytokine,
                             aliases)) %>%
    select(-cytokine, cytokine=aliases)
active_cytokines

# Filter liana according to the positve cytokine activities alone
liana_filt <- liana_res %>%
    left_join(active_cytokines, by=c("ligand"="cytokine", "target"="celltype")) %>%
    filter(active_flag == 1)

# plot
liana_dotplot(liana_filt,
              source_groups = c("B"),
              target_groups = c("NK", "CD8 T"))
# Note if a dot is missing, it simply means that the interaction is filtered out
# i.e. the ligand or receptor are not expressed in at least 10% of the cells
