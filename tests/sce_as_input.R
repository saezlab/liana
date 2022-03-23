# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))
sce <- Seurat::as.SingleCellExperiment(seurat_object)
SingleCellExperiment::colLabels(sce) <- Idents(seurat_object)


op_resource = select_resource("OmniPath")[[1]] %>%
    decomplexify()

assay = "RNA"

#
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
entity_genes = union(transmitters$gene,
                     receivers$gene)


entity_genes


# Filter to LigRec and scale
seurat_object <- seurat_object[rownames(seurat_object) %in% entity_genes]
seurat_object <- Seurat::ScaleData(seurat_object, features = entity_genes)

sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)
colLabels(sce) <- SeuratObject::Idents(seurat_object)
sce@assays@data$scaledata <- seurat_object@assays[[assay]]@scale.data


res1 <- liana_pipe(sce = liana_prep(sce = seurat_object),
                   op_resource = select_resource("OmniPath")[[1]] %>%
                       decomplexify())

liana_prep(sce)
liana_prep(seurat_object)


sce <- sce[rowSums(counts(sce)) > 0,
           colSums(counts(sce)) > 0]



astro <- readRDS("~/Downloads/astro.data_cortex.rds")
astro <- astro %>% Seurat::NormalizeData()

xd <- liana_wrap(astro,
                 method = "sca",
                 resource = "custom",
                 external_resource = generate_orthologs(op_resource,
                                                        symbols_dict),
                 idents_col = "Type"
           )



rownames(astro@assays$RNA@counts) <- toupper(rownames(astro))
rownames(astro@assays$RNA@data) <- toupper(rownames(astro))


xd2 <- liana_wrap(astro, idents_col = "Type", method="sca")




require(GENIE3)
GENIE3::GENIE3(exprMatrix = as.matrix(astro@assays$RNA@counts),
               regulators = rownames(astro@assays$RNA@counts)
               )

