liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))


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

# liana Pipe Output ----
pipe_out <- liana_pipe(seurat_object,
                       op_resource = select_resource("OmniPath")[[1]] %>%
                           liana:::decomplexify())



pipe_out %>%
    select(source, target,
           ligand, receptor,
           ligand.expr, receptor.expr) %>%
    mutate(ligand.soft = softmax(ligand.expr)) %>%
    mutate(receptor.soft = softmax(receptor.expr)) %>%
    mutate(comm_score = pmap_dbl(.l = list(.$ligand.soft, .$receptor.soft),
                             .f = function(l, r){
                                 m <- matrix(cbind(l, r), nrow = 1)
                                 norm(m, type = "2")
                                 }
                             )) %>%
    arrange(desc(comm_score))









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
cellcall_out <- readRDS("~/Downloads/cellcall_out.RDS")
gene_means <- cellcall_out@data$expr_mean %>% as_tibble(rownames="gene")

# Start of Function
my_Expr <- cellcall_out@data$withoutlog
colnames(my_Expr) <- as.character(cellcall_out@meta.data$celltype)
my_Expr[1:4,1:4]
detect_gene <- rownames(my_Expr)



#
object = cellcall_out
pValueCor=0.05
CorValue=0.1
topTargetCor=1
method="weighted"
p.adjust=0.05
use.type="median"
probs = 0.9
Org = 'Homo sapiens'
IS_core = TRUE

# Read LR-TF-Pathway df
f.tmp <- system.file("extdata", "new_ligand_receptor_TFs.txt", package="cellcall")
triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcall")
triple_relation
triple_relation$pathway_ID <- NULL
print(triple_relation[1:4,])

# TF-regulons
target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
tmp_complex_symbol <- triple_relation$Receptor_Symbol[grep(",",triple_relation$Receptor_Symbol)] %>% unique() %>% str_split(",") %>% unlist %>% unique()

all.gene.needed <- unique(as.character(c(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, target_relation$TF_Symbol, target_relation$Target_Symbol,tmp_complex_symbol)))
# triple_relation[1:4,1:4]
# target_relation[1:4,1:4]
my_Expr <- object@data$withoutlog
colnames(my_Expr) <- as.character(object@meta.data$celltype)
my_Expr[1:4,1:4]



expr_set <- my_Expr[intersect(detect_gene, all.gene.needed),]
detect_gene <- rownames(expr_set)
cell_type = unique(colnames(expr_set))
expr.fc <- object@data$withoutlog[detect_gene,]
colnames(expr.fc) <- colnames(expr_set)

rm(list=c("object"))

complex_matrix <- matrix(ncol = length(colnames(expr_set)))
complex_matrix <- as.data.frame(complex_matrix)
colnames(complex_matrix) <- colnames(expr_set)
myrownames <- c()

expr_set <- expr_set[apply(expr_set, 1, function(x){sum(x!=0)})>0,]
detect_gene <- rownames(expr_set)
# expr_set[1:4,1:4]

print("step1: compute means of gene")
expr_mean <- matrix(nrow = nrow(expr_set), ncol = length(cell_type))
myColnames <- c()
for (i in 1:length(cell_type)) {
    myCell <- cell_type[i]
    myMatrix <- expr_set[,colnames(expr_set)==myCell,drop=F]
    if(use.type=="mean"){
        myMatrix_mean <- as.numeric(apply(myMatrix, 1, mean))
    }else if(use.type=="median"){
        quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
            quantile(x, probs = probs,names=FALSE)
        }))
        mean.tmp <- rowMeans(myMatrix)
        mean.tmp[which(quantil.tmp==0)]<-0
        myMatrix_mean <- mean.tmp
    }
    expr_mean[,i] <- myMatrix_mean
    myColnames <- c(myColnames, myCell)
    # print(myCell)
}
expr_mean <- data.frame(expr_mean)
colnames(expr_mean) <- myColnames
rownames(expr_mean) <- rownames(expr_set)

expr_mean <- expr_mean[apply(expr_mean, 1, function(x){sum(x!=0)})>0,]
detect_gene <- rownames(expr_mean)


#
ligand_symbol <- unique(triple_relation$Ligand_Symbol)
softmax_ligand <- expr_mean[intersect(ligand_symbol, detect_gene),]
colnames(softmax_ligand) <- colnames(expr_mean)
rowCounts <- rowSums(softmax_ligand)

softmax_ligand <- do.call(rbind,lapply(1:nrow(softmax_ligand), function(i){
    softmax_ligand[i,]/rowCounts[i]
}))


# Notes:
# CellCall does not use the softmax function :D, it simply divides by the rowSums
# i.e. each ligand is divided by the max of the same ligand.



# Create seurat - TEST IN LIANA ----------
seurat_test <- CreateSeuratObject(CreateAssayObject(counts=cellcall_out@data$count)) %>%
    AddMetaData(cellcall_out@meta.data)
Idents(seurat_test) <- seurat_test@meta.data$celltype

#
liana_res <- liana_pipe(seurat_test,
                        op_resource = select_resource("OmniPath")[[1]] %>%
                            decomplexify())

liana_res %<>%
    group_by(ligand) %>%
    mutate(ligand.max = max(ligand.sum)) %>%
    ungroup() %>%
    mutate(ligand.softmax = ligand.expr/ligand.max) %>%
    group_by(receptor) %>%
    mutate(receptor.max = max(receptor.sum)) %>%
    ungroup() %>%
    mutate(receptor.softmax = receptor.expr/receptor.max) %>%
    rowwise() %>%
    mutate(cellcall.score = mean(c(ligand.softmax,receptor.softmax)))



cellcall_lr <- cellcall_out@data$expr_l_r %>% as_tibble(rownames="lr")


#### Resource test
cellcall_lr <- read.table("/home/dbdimitrov/R/x86_64-pc-linux-gnu-library/4.0/cellcall/extdata/new_ligand_receptor_TFs.txt", header=TRUE)
# * extended table is also available but the number does not match with what's reported in the paper

cellcall_lr %>%
    select(Ligand_Symbol, Receptor_Symbol) %>%
    distinct() %>%
    nrow()

