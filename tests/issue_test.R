seurat_object <- readRDS("~/Repos/ligrec_decouple/data/input/spatial/fishes/seqFISH_seurat.RDS")

liana_res <- liana_wrap(seurat_object, method="call_sca", assay="SCT")



op_resource <- readRDS("~/Repos/ligrec_decouple/data/input/murine_omnipath.RDS")
.format = TRUE
assay = "RNA"
assay.type = "logcounts"

if(assay.type=="logcounts"){
    assay.type = "data"
}

# Format OmnipathR resource
if(!is.null(op_resource)){
    op_resource %<>% sca_formatDB
} else{
    if(file.exists(system.file(package = "liana", "LRdb.rda"))){
        load(system.file(package = "liana", "LRdb.rda"))
        op_resource <- LRdb
    } else{
        stop("Could not locate LRdb.rda")
    }
}

# Prepare data from Seurat object
input_data <-
    Seurat::GetAssayData(seurat_object,
                         assay = assay,
                         slot = assay.type)
labels <- Seurat::Idents(seurat_object)

# Compute interactions between cell clusters
signal <- SCAomni::cell_signaling(data = input_data,
                                  genes = row.names(input_data),
                                  cluster = as.numeric(labels),
                                  c.names = levels(Idents(seurat_object)),
                                  species = 'homo sapiens',
                                  LRdb = op_resource,
                                  int.type="autocrine", # includes both para and auto...
                                  write = FALSE,
                                  verbose = FALSE#,#
                                  # ...
)


# Compute intercellular gene networks
sca_res <- SCAomni::inter_network(data = input_data,
                                  signal = signal,
                                  genes = row.names(input_data),
                                  cluster = as.numeric(labels),
                                  c.names = levels(Idents(seurat_object)),
                                  write = FALSE
)
if (.format) {
    sca_res %<>% FormatSCA
}
