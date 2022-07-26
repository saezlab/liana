library(ExperimentHub)
eh <- ExperimentHub()
query(eh, "Kang")

# Get Data
(sce <- eh[["EH2259"]])

sce
dim(sce)

# Preprocessing
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)
# names are not ok?
rownames(sce) <- gsub(rownames(sce), pattern = "_ENS*", replacement = "")

# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)
qc


# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)


# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# compute sum-factors & normalize
# Finally, we use logNormCounts to calculate log2-transformed normalized
# expression values by dividing each count by its size factor,
# adding a pseudo-count of 1, and log-transforming
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# Note that, in this workflow, expression values are used for visualization only,
# and that differential analyses are performed on pseudobulks (section 3.3)
# or the count data directly (section 3.4).

# # Alternatively, expression values could be obtained via vst
# # (variance stabilizing transformation) from the sctransform package (Hafemeister and Satija 2019),
# # which returns Pearson residuals from a regularized negative binomial regression
# # model that can be interpreted as normalized expression values:
#
# library(sctransform)
# assays(sce)$vstresiduals <- vst(counts(sce), verbosity = FALSE)$y


# Muscat Data preparation ----
require(tidyverse)
library(muscat)

sce$id <- paste0(sce$stim, sce$ind)
# Prepare SCE for DS analysis
(sce <- prepSCE(sce,
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

# # compute UMAP using 1st 20 PCs
# sce <- runUMAP(sce, pca = 20)
#
# # wrapper to prettify reduced dimension plots
# .plot_dr <- function(sce, dr, col){
#     plotReducedDim(sce, dimred = dr, colour_by = col) +
#         guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
#         theme_minimal() + theme(aspect.ratio = 1)
# }
#
# # downsample to max. 100 cells per cluster
# cs_by_k <- split(colnames(sce), sce$cluster_id)
# cs100 <- unlist(sapply(cs_by_k, function(u)
#     sample(u, min(length(u), 100))))
#
# # plot t-SNE & UMAP colored by cluster & group ID
# for (dr in c("TSNE", "UMAP"))
#     for (col in c("cluster_id", "group_id"))
#         plot(.plot_dr(sce[, cs100], dr, col))

# Pseudobulk
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

# # prop
# pb.prop <- aggregateData(sce,
#                          assay = "counts", fun = "prop.detected",
#                          by = c("cluster_id", "sample_id"))
#
# pb.prop@assays@data$`B cells`


#  Pseudobulk-level MDS plot
pb_mds <- pbMDS(pb) +
    scale_shape_manual(values = c(17, 4)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# change point size & alpha
pb_mds$layers[[1]]$aes_params$size <- 5
pb_mds$layers[[1]]$aes_params$alpha <- 0.6
plot(pb_mds)


# Sample-level analysis: Pseudobulk methods ----

# What happens inside ====
# Define contrast
contrasts <- list(
    "stim" = "group_idstim - group_idctrl"
)
design = model.matrix(~ 0 + group_id, data=as.data.frame(colData(pb)))

# Contrast matrix
cont_matrix <- limma::makeContrasts(contrasts = contrasts,
                                    levels = design)

# Check min cells
muscat:::.n_cells(pb)


# Filter by min cells and do limma
# this function is modified for now should install muscat from saezlab
t_res <- muscat::pbDS(pb,
                      method = "limma-voom",
                      design = design,
                      contrast = cont_matrix,
                      filter = "both",
                      min_cells = 10,
                      # edgeR args
                      min.prop = 0.5,
                      min.count = 5,
                      min.total.count = 10)


# access results table for 1st comparison
tbl <- t_res$table[[1]]
# one data.frame per cluster
names(tbl)

# Get all Possible Cluster pair combinations ---
pairs <- expand_grid(source = names(tbl),
                     target = names(tbl))

op_resource <- liana::select_resource("Consensus")[[1]] %>%
    liana::decomplexify() %>%
    select(source, target,
           source_genesymbol, target_genesymbol,
           source_genesymbol_complex, target_genesymbol_complex)

# Resource Format
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value) %>%
    distinct()
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value) %>%
    distinct()

# Find Markers and Format (Example for join)
cols <- c("gene", "logFC", "AveExpr", "stat",
          "p_val", "p_adj.loc", "p_adj.glb")


cluster_markers <- t_res$table[[1]] %>%
    map(function(cluster){
        cluster %>%
            dplyr::rename(stat=t) %>% # !
            as_tibble() %>%
            select(all_of(cols))
    })


lr_res <- pairs %>%
    pmap(function(source, target){
        source_stats <- ligrec_degformat(entity = transmitters,
                                         ligand_receptor = "ligand",
                                         columns = cols,
                                         cluster_markers = cluster_markers[[source]])
        target_stats <- ligrec_degformat(entity = receivers,
                                         ligand_receptor = "receptor",
                                         columns = cols,
                                         cluster_markers = cluster_markers[[target]])

        op_resource %>%
            select(ligand = source_genesymbol,
                   receptor = target_genesymbol) %>%
            left_join(source_stats, by = "ligand") %>%
            left_join(target_stats, by = "receptor") %>%
            distinct() %>%
            mutate(source = source,
                   target = target) %>%
            na.omit()
    }) %>%
    bind_rows()


# Join complexes (recomplexify) to lr_res
cmplx <- op_resource %>%
    select(
        ligand = source_genesymbol,
        ligand.complex = source_genesymbol_complex,
        receptor = target_genesymbol,
        receptor.complex = target_genesymbol_complex
    )

# Recomplexify
lr_res %<>%
    left_join(., cmplx,
              by=c("ligand", "receptor")) %>%
    distinct()
lr_res %<>%
    recomplexify(
        columns = c("ligand.logFC", "ligand.AveExpr",
                    "receptor.logFC", "receptor.AveExpr",
                    "ligand.p_val", "receptor.p_val",
                    "ligand.p_adj.loc", "receptor.p_adj.loc",
                    "ligand.p_adj.glb", "receptor.p_adj.glb",
                    "ligand.stat", "receptor.stat"),
        complex_policy="mean0") %>%
    ungroup() %>%
    filter(ligand.AveExpr!=0 & receptor.AveExpr!=0)

# the golden mean
golden_mean <- function(vec){ mean(vec) / 1.618 ^ sd(vec) }

# Calculate mean stat
lr_res %<>%
    rowwise() %>%
    mutate(lr_stat = mean(c(ligand.stat, receptor.stat))) %>%
    mutate(lr_std = golden_mean(c(ligand.stat, receptor.stat))) %>%
    dplyr::rename_with(~gsub(x=.x,
                             pattern="p_val",
                             replacement = "pval"),
                       ends_with("p_val")) %>%
    dplyr::rename_with(~gsub(x=.x,
                             pattern="AveExpr",
                             replacement = "expr"),
                       ends_with("AveExpr")) %>%
    dplyr::select(source, target, ligand, receptor,
                  lr_stat, lr_std,
                  everything()) %>%
    arrange(desc(abs(lr_stat)))




