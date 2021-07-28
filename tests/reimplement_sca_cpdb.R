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



# LRscore ----
# SCA returns only paracrine interactions...
# albeit the default is int.type=c("paracrine", "autocrine"), it works exclusively
# for some unknown reason, also some interactions simply disappear and it's impossible to track where?
global_mean <- test_sce@assays@data$logcounts %>%
    .[Matrix::rowSums(.)>0,]
global_mean <- sum(global_mean)/(nrow(global_mean)*ncol(global_mean))


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
    mutate(LRscore = LRscore(ligand.expr, receptor.expr, global_mean)) %>%
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



# CPDB ----
lr_cpdb <- lr_res %>%
    filter(receptor.prop >= 0.1 & ligand.prop >= 0.1) %>%
    select(-ends_with(c(".pval", "scaled", ".FDR", "stat"))) %>%
    rowwise() %>%
    mutate(lr.mean = mean(ligand.expr, receptor.expr, trim=0))


# taken from CellChat
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
    percent <- Matrix::nnzero(x)/length(x) # not sure how useful?
    if (percent < trim) {
        return(0)
    } else {
        return(mean(x, na.rm = na.rm, trim=trim))
    }
}

# columns
group <- colLabels(test_sce)

trunc_mean <- aggregate(t(as.matrix(test_sce@assays@data$counts)),
                        list(group),
                        FUN=thresholdedMean, trim=0.05) %>%
    as_tibble() %>%
    rename(celltype = Group.1) %>%
    pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from=celltype, id_cols=gene,values_from=value) %>%
    column_to_rownames("gene")




require(forcats)
fct_shuffle(group)




# shuffle columns
clusts <- colLabels(test_sce) %>%
    as_tibble(rownames = "cell")




shuffled_clusts <-
    map(1:1000, function(perm){
        clusts %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })


lrs <- lr_cpdb %>%
    select(ligand,receptor,source,target)

# parallelize
require(furrr)
st <- Sys.time()
plan(multisession, workers = 8)
ff <- furrr::future_map(shuffled_clusts, function(clust){
    aggregate(t(as.matrix(test_sce@assays@data$counts)),
              list(clust),
              FUN=thresholdedMean,
              trim=0.05) %>%
        as_tibble() %>%
        rename(celltype = Group.1) %>%
        pivot_longer(-celltype, names_to = "gene") %>%
        tidyr::pivot_wider(names_from=celltype,
                           id_cols=gene,
                           values_from=value) %>%
        column_to_rownames("gene")
    })
Sys.time() - st




st <- Sys.time()
pp <- map(shuffled_clusts, function(clust){
    aggregate(t(as.matrix(test_sce@assays@data$counts)),
              list(clust), FUN=thresholdedMean, trim=0.05) %>%
        as_tibble() %>%
        rename(celltype = Group.1) %>%
        pivot_longer(-celltype, names_to = "gene") %>%
        tidyr::pivot_wider(names_from=celltype, id_cols=gene,values_from=value) %>%
        column_to_rownames("gene")
    })
st - Sys.time()




rand_cpdb <-
    ff %>% map(function(perm_means){
        lrs %>%
            join_means(means = perm_means,
                       source_target = "source",
                       entity = "ligand",
                       type = "trunc") %>%
            join_means(means = perm_means,
                       source_target = "target",
                       entity = "receptor",
                       type = "trunc") %>%
            mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
    })

to_check <- lr_cpdb %>%
    select(ligand, receptor, source, target, non_random=lr.mean)

bind_perms <- rand_cpdb %>%
    bind_rows() %>%
    left_join(to_check, by = c("ligand", "receptor", "source", "target"))


pvals_df <- bind_perms %>%
    group_by(ligand, receptor, source, target) %>%
    summarise(pval = 1 - sum(non_random >= lr_mean)/1000)





# join trunc mean
lr_cpdb2 <- lr_cpdb %>%
    join_means(means = trunc_mean,
               source_target = "source",
               entity = "ligand",
               type = "trunc") %>%
    join_means(means = trunc_mean,
               source_target = "target",
               entity = "receptor",
               type = "trunc") %>%
    mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))



lr_cp <- lr_cpdb2 %>%
    select(ligand, receptor, source, target, lr.mean) %>%
    left_join(pvals_df)







res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key=NULL,
                     n_perms=1000,
                     threshold=0.01,
                     seed=as.integer(1004))


