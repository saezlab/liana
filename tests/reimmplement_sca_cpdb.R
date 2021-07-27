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
test_sce@assays@data$scaledata <- seurat_object@assays$RNA@scale.data
colLabels(test_sce) <- Seurat::Idents(seurat_object)


xx <- test_sce@assays@data$logcounts
xx <- xx[Matrix::rowSums(xx)>0,]
med <- sum(xx)/(nrow(xx)*ncol(xx))
med


global_mean <- test_sce@assays@data$logcounts %>%
    .[Matrix::rowSums(.)>0,]
global_mean <- sum(global_mean)/(nrow(global_mean)*ncol(global_mean))



# LRscore
# SCA returns only paracrine interactions...
# albeit the default is int.type=c("paracrine", "autocrine"), it works exclusively
# for some unknown reason, also some interactions simply disappear and it's impossible to track where?
exp1 <- call_sca(op_resource = select_resource("OmniPath")[[1]],
                 seurat_object = seurat_object,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(0),
                 int.type = c("paracrine", "autocrine"))

res1 <- lr_res %>%
    select(-ends_with("scaled")) %>%
    select(-ends_with(".pval")) %>%
    mutate(LRscore = LRscore(ligand.expr, receptor.expr, med)) %>%
    arrange(desc(LRscore))



xd <- get_sca(op_resource = select_resource("OmniPath")[[1]],
              seurat_object = seurat_object)


res1 <- liana_wrap(seurat_object,
                   method = c('sca','call_sca'),
                   resource = c('OmniPath'),
                   sca.params=list(logFC=0,
                                   s.score=0,
                                   int.type="autocrine",
                                   tol=0))

res1
