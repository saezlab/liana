# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
op_resource <- select_resource("OmniPath")[[1]]


lr_res <- liana_pipe(seurat_object,
                     op_resource)
lr_res


# OP format
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)

entity_genes <- union(transmitters$gene, receivers$gene)


# Convert to SCE
seurat_object <- seurat_object[rownames(seurat_object) %in% entity_genes]
seurat_object <- Seurat::ScaleData(seurat_object, features = entity_genes)
test_sce <- Seurat::as.SingleCellExperiment(seurat_object,
                                            assay="RNA")
test_sce@assays@data$logcounts <- seurat_object@assays$RNA@scale.data
colLabels(test_sce) <- Seurat::Idents(seurat_object)


# Connectome ----
conn_exp <- readRDS(file.path(liana_path, "testdata",
                              "output", "conn_res.RDS"))

conn_res <- lr_res %>%
    rowwise() %>%
    mutate(conn_score = mean(c(ligand.scaled, receptor.scaled))) %>%
    select(source, ligand, ligand.scaled,
           target, receptor, receptor.scaled, conn_score)

scaled <- scuttle::summarizeAssayByGroup(test_sce,
                                         ids = colLabels(test_sce),
                                         assay.type = "logcounts")
scaled <- scaled@assays@data$mean  %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble()

scaled_s <- Seurat::AverageExpression(seurat_object, slot = "scale.data") %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble()



lr_res




# NATMI ----

# expected
natmi_exp <- readRDS(file.path(liana_path, "testdata",
                               "output", "natmi_res.RDS"))

natmi_res <- lr_res %>%
    select(source, ligand, ligand.count, ligand.sum,
           target, receptor, receptor.count, receptor.sum) %>%
    rowwise() %>%
    mutate(natmi_score = ((ligand.count*(ligand.sum^-1))) *
               ((receptor.count*(receptor.sum^-1)))) %>%
    arrange(desc(natmi_score)) %>%
    filter(natmi_score!=0) %>%
    arrange(ligand, receptor, source, target) %>%
    unite(ligand, receptor, source, target,
          col = "interaction", sep = "_")

natmi_exp <- natmi_exp %>%
    arrange(ligand, receptor, source, target) %>%
    unite(ligand, receptor, source, target,
          col = "interaction", sep = "_")

not_there <- natmi_exp %>%
    filter(!(interaction %in% natmi_res$interaction))




# Correlation + p-value ----
corr_pairs <- correlatePairs(test_sce,
                             subset.row =
                                 which(rownames(test_sce) %in% entity_genes)
                             ) %>%
    as_tibble()

corr_pairs <- corr_pairs %>%
    select(gene1=gene2,
           gene2=gene1,
           everything()) %>%
    bind_rows(corr_pairs) %>%
    select(gene1,
           gene2,
           rho,
           corr.FDR=FDR)


corr_score <- lr_res %>%
    left_join(
        corr_pairs,
        by=c("ligand"="gene1",
             "receptor"="gene2")
        ) %>%
    distinct() %>%
    filter(ligand.FDR <= 0.05,
           receptor.FDR <= 0.05) %>%
    arrange(desc(rho))


# logFC product
markers <- FindAllMarkers(seurat_object)

markers_sce <- scran::findMarkers(test_sce,
                                  groups = colLabels(test_sce),
                                  direction = "any",
                                  full.stats = TRUE,
                                  test.type = "t",
                                  pval.type = "all",
                                  assay.type = "counts")

markers_sce <- scran::pairwiseTTests(test_sce,
                                  groups = colLabels(test_sce),
                                  direction = "any",
                                  assay.type = "counts")
markers_sce

b_nk_folds <- Seurat::FoldChange(seurat_object, ident.1="B", ident.2 ="NK")

# Get Avg Per Cluster (data assay)
means <- scuttle::summarizeAssayByGroup(test_sce,
                                        ids = colLabels(test_sce),
                                        assay.type = "logcounts")
means <- means@assays@data$mean

