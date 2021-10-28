softmax <- function(par){
    n.par <- length(par)
    par1 <- sort(par, decreasing = TRUE)
    Lk <- par1[1]
    for (k in 1:(n.par-1)) {
        Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk)))
    }
    val <- exp(par - Lk)
    return(val)
}

# Example 1
vec <- c(-1,2,1,-3)
sm <- softmax(vec)
print(sm)


z <- c(1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0)
softmax <- exp(z)/sum(exp(z))
softmax


# Test CellCall
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
require(cellcall)
f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellcall")
load(f.tmp)

#
cc_object <- cellcall::CreateObject_fromSeurat(seurat_object)

mt <- TransCommuProfile(object = cc_object,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')



## gene expression stored in the variable in.content
dim(in.content)
in.content[1:4, 1:4]
table(str_split(colnames(in.content), "_", simplify = T)[,2])


mt <- CreateNichConObject(data=in.content, min.feature = 3,
                          names.field = 2,
                          names.delim = "_",
                          source = "TPM",
                          scale.factor = 10^6,
                          Org = "Homo sapiens",
                          project = "Microenvironment")

mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')
saveRDS(mt, "~/Downloads/cellcall_out.RDS")

# from CellCall object to Seurat
seurat_test <- CreateSeuratObject(CreateAssayObject(counts=mt@data$count))
seurat_test <- AddMetaData(seurat_test, mt@meta.data)
Idents(seurat_test) <- seurat_test@meta.data$celltype

#
liana_res <- liana_pipe(seurat_test,
                        op_resource = select_resource("OmniPath")[[1]] %>%
                            decomplexify())
liana_res %<>%
    mutate(ligand.softmax = (ligand.expr)/ligand.sum)

liana_res %>%
    filter(source=="ST")

cellcall_ligand <- mt@data$softmax_ligand
cellcall_ligand %>% head()
