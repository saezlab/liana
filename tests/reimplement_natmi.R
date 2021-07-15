# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
op_resource <- select_resource("OmniPath")[[1]]

# expected
natmi_exp <- readRDS(file.path(liana_path, "testdata",
                               "output", "natmi_res.RDS"))


require(SingleCellExperiment)
require(scuttle)
require(scran)

lr_res <- liana_pipe(seurat_object,
                     op_resource,
                     test.type = "t")


# NATMI
natmi_res <- lr_res %>%
    select(source, ligand, ligand.avg, ligand.sum,
           target, receptor, receptor.avg, receptor.sum) %>%
    rowwise() %>%
    mutate(natmi_specificity = ((ligand.avg*(ligand.sum^-1))) *
               ((receptor.avg*(receptor.sum^-1)))) %>%
    arrange(desc(natmi_specificity)) %>%
    filter(natmi_specificity!=0) %>%
    arrange(ligand, receptor, source, target) %>%
    unite(ligand, receptor, source, target,
          col = "interaction", sep = "_")


# t-value product
tvalue_lr <- lr_res %>%
    rowwise() %>%
    filter(ligand.FDR <= 0.05 && receptor.FDR <= 0.05) %>%
    mutate(t_score = ligand.stat * receptor.stat)





natmi_exp <- natmi_exp %>%
    arrange(ligand, receptor, source, target) %>%
    unite(ligand, receptor, source, target,
          col = "interaction", sep = "_")


not_there <- natmi_exp %>%
    filter(!(interaction %in% natmi_res$interaction))




# Convert to SCE
test_sce <- Seurat::as.SingleCellExperiment(seurat_object,
                                            assay="RNA")
colLabels(test_sce) <- Seurat::Idents(seurat_object)
# test_sce <- scuttle::logNormCounts(test_sce)



# Get Avg Per Cluster (data assay)
means <- scuttle::summarizeAssayByGroup(test_sce,
                                        ids = colLabels(test_sce),
                                        assay.type = "counts")
means <- means@assays@data$mean

natmi_res <- liana_pipe(seurat_object,
                        op_resource,
                        test.type = "t") %>%
    select(source, ligand, ligand.expr,
           target, receptor, receptor.expr)


# Get Average per Gene (i.e. global average)
means_avg <- means@assays@data$mean %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    rowwise("gene") %>%
    mutate(mean.sum = mean(c_across(where(is.numeric)))) %>%
    select(gene, mean.sum) %>%
    dplyr::rename(ligand = gene,
                  ligand.sum = mean.sum) %>%
    distinct()




