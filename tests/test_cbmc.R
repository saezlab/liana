liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

Idents(seurat_object) %<>%
    as.numeric() %>% as.character() %>% as.factor()

seurat_object@meta.data$seurat_annotations <- Idents(seurat_object)
seurat_object@meta.data

liana_wrap(seurat_object,
           method="call_natmi")


# Load and and format Data
seurat_object <- readRDS("~/Repos/ligrec_decouple/data/input/cbmc_seurat.rds")

liana_wrap(seurat_object,
           method = "connectome")


# Call liana
sca_results <-
    call_sca(op_resource = NULL,
             seurat_object,
             assay = 'RNA',
             .format = TRUE,
             s.score = 0,
             logFC = log2(1.5))

italk_results <- call_italk(op_resource = NULL,
                            seurat_object,
                            assay = 'RNA',
                            .format = TRUE,
                            .DE = TRUE)

conn_res <- call_connectome(seurat_object = seurat_object,
                            op_resource = select_resource("OmniPath")[[1]],
                            min.cells.per.ident = 1,
                            p.values = TRUE,
                            calculate.DOR = FALSE,
                            .format = TRUE,
                            assay = 'RNA')


cmbc_cc <- liana_wrap(
    seurat_object = seurat_object,
    resource = c("OmniPath"),
    method = "connectome"
)

op_resource <- select_resource("OmniPath")[[1]]

lr_db <- conn_formatDB(op_resource)

# scale genes to ligands and receptors available in the resource
connectome.genes <- union(lr_db$source_genesymbol, lr_db$target_genesymbol)

genes <- connectome.genes[connectome.genes %in% rownames(seurat_object)]
seurat_objectx <- Seurat::ScaleData(seurat_object,
                                    features = genes,
                                    nbrOfWorkers = 1)

seurat_object@assays$RNA@scale.data
seurat_objectx@assays$RNA@scale.data
