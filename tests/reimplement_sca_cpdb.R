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


# trunc mean
trunc_mean <- aggregate(t(as.matrix(test_sce@assays@data$counts)),
                        list(colLabels(test_sce)),
                        FUN=mean, trim=0.00) %>%
    as_tibble() %>%
    rename(celltype = Group.1) %>%
    pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from=celltype, id_cols=gene,values_from=value) %>%
    column_to_rownames("gene")


# filter by proportion and join trunc mean
lr_cpdb <- lr_res %>%
    filter(receptor.prop >= 0.01 & ligand.prop >= 0.01) %>% # to be moved
    select(-ends_with(c(".pval", "scaled", ".FDR", "stat"))) %>%
    join_means(means = trunc_mean,
               source_target = "source",
               entity = "ligand",
               type = "trunc") %>%
    join_means(means = trunc_mean,
               source_target = "target",
               entity = "receptor",
               type = "trunc") %>%
    rowwise() %>%
    mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))



st <- Sys.time()
cpdb_res <- cpdb_score(lr_res = lr_cpdb,
                       sce_mat = t(as.matrix(test_sce@assays@data$counts)),
                       nperms = 10000,
                       seed = 1234,
                       trim = 0.05,
                       parallelize=TRUE,
                       workers = 4)
Sys.time() - st

st <- Sys.time()
res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key=NULL,
                     n_perms=10000,
                     threshold=0.05,
                     seed=as.integer(1004))
Sys.time() - st




## Generate LR distributions
# save expression matrix
sce_matrix <- t(as.matrix(test_sce@assays@data$counts))

nperms <- 1000
seed <- 1234

# shuffle columns
set.seed(seed)
shuffled_clusts <-
    map(1:nperms, function(perm){
        colLabels(test_sce) %>%
            as_tibble(rownames = "cell") %>%
            slice_sample(prop=1, replace = FALSE) %>%
            deframe()
    })


# keep all lrs
lr_og <- lr_cpdb %>%
    select(ligand, receptor, source, target, og_mean = lr_mean)




# future_map
require(future)
require(furrr)
future::plan(multisession, workers = 4)
perm <- furrr::future_map(shuffled_clusts,
                          ~liana:::cpdb_permute(
                              col_labels = .x,
                              sce_matrix=sce_matrix,
                              trim = 0,
                              lr_res)
                          )

pvals_df <- perm %>%
    bind_rows() %>%
    left_join(to_check, by = c("ligand", "receptor", "source", "target")) %>%
    group_by(ligand, receptor, source, target) %>%
    mutate(lr_mean = na_if(lr_mean, 0)) %>%
    summarise(pval = 1 - (sum(og_mean >= lr_mean)/nperms)) #nperm


lr_cp <- lr_cpdb %>%
    select(ligand, receptor, source, target, lr_mean) %>%
    left_join(pvals_df, by = c("ligand", "receptor", "source", "target"))



res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key=NULL,
                     n_perms=1000,
                     threshold=0.01,
                     seed=as.integer(1004))



















xx <- all_lrs %>%
    join_means(means = pp[[1]],
               source_target = "source",
               entity = "ligand",
               type = "trunc") %>%
    join_means(means = pp[[1]],
               source_target = "target",
               entity = "receptor",
               type = "trunc") %>%
    rowwise() %>%
    mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))





# parallelize
require(furrr)
st <- Sys.time()
plan(multisession, workers = 4)
ff <- furrr::future_map(shuffled_clusts, function(clust){
    aggregate(t(as.matrix(test_sce@assays@data$counts)),
              list(clust), FUN=mean, trim=0) %>%
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
              list(clust), FUN=mean, trim=0.05) %>%
        as_tibble() %>%
        rename(celltype = Group.1) %>%
        pivot_longer(-celltype, names_to = "gene") %>%
        tidyr::pivot_wider(names_from=celltype, id_cols=gene,values_from=value) %>%
        column_to_rownames("gene")
    })
st - Sys.time()




rand_cpdb <-
    ff %>% map(function(perm_means){
        all_lrs %>%
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
    select(ligand, receptor, source, target, non_random=lr_mean)

bind_perms <- rand_cpdb %>%
    bind_rows() %>%
    left_join(to_check, by = c("ligand", "receptor", "source", "target"))


pvals_df <- bind_perms %>%
    group_by(ligand, receptor, source, target) %>%
    summarise(pval = 1 - sum(non_random >= lr_mean)/1000)








lr_cp <- lr_cpdb2 %>%
    select(ligand, receptor, source, target, lr.mean) %>%
    left_join(pvals_df, by = c("ligand", "receptor", "source", "target"))




res1 <- call_squidpy(seurat_object = seurat_object,
                     op_resource = select_resource("OmniPath"),
                     cluster_key=NULL,
                     n_perms=1000,
                     threshold=0.01,
                     seed=as.integer(1004))


