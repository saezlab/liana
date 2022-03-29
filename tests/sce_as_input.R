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


# Test log and antilog
# input
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                                   "input", "testdata.rds"))

test_external=FALSE # Set to true in any R/ script to test external pipes

# convert and save the Seurat object
seurat_conv <- liana_prep(seurat_object)[1:50,]

normmat <- scuttle::normalizeCounts(seurat_conv, log=TRUE)
libmat  <- scuttle::normalizeCounts(seurat_conv, log=FALSE)
logmat <- log2(counts(seurat_conv)@x + 1)
rawmat <- counts(seurat_conv)


normmat@x[1:10]
libmat@x[1:10]
rawmat@x[1:10]

# Should be equal to libmat

antilog1m(slot(normmat, "x"))
normmat@x[1:10]

normmat

normmat_og <- normmat

assay.type = "logcounts"
sce = seurat_conv
antilogged <- antilog1m(slot(exec(assay.type, sce), "x"))
sce@assays@data[[assay.type]]@x <- antilogged

get_log2FC(sce, "logcounts")




round(antilog1m(normmat@x)[1:10], 5) == round(libmat@x[1:10], 5)

.antilog1m <- function(x, base=2){base ^ (x) - 1}

.get_normcounts <- function(sce, assay.type){
    antilogged <- .antilog1m(slot(exec(assay.type, sce), "x"))

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

.get_normcounts(sce, "logcounts")

xx <- rbinom(100, 100, 0.1)
xx

antilog1m(log2(xx + 1))
