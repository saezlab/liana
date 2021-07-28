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
    filter(receptor.prop >= 0.01 & ligand.prop >= 0.1) %>%
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

trunc_mean <- aggregate(t(as.matrix(test_sce@assays@data$counts)), list(group), FUN=thresholdedMean, trim=0.05)
trunc_mean <- t(trunc_mean[,-1])
colnames(trunc_mean) <- levels(group)
head(trunc_mean)

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

require(furrr)
st <- Sys.time()
plan(multisession, workers = 8)
ff <- furrr::future_map(shuffled_clusts, function(clust){
    aggregate(t(as.matrix(test_sce@assays@data$counts)),
              list(clust),
              FUN=thresholdedMean,
              trim=0.05)
    })
Sys.time() - st


st <- Sys.time()
pp <- map(shuffled_clusts, function(clust){
    aggregate(t(as.matrix(test_sce@assays@data$counts)), list(clust), FUN=thresholdedMean, trim=0.05)
    })
st - Sys.time()




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
    mutate(lr_mean = mean(ligand.trunc, receptor.trunc))


permutation <- replicate(100, sample.int(90, size = 90))

groupboot <- group[permutation[, 1]]
groupboot2 <- group[permutation[, 2]]



require(tidymodels)
tidymodels_prefer()


# resample1 <- lr_cpdb2 %>% ungroup() %>% permutations(lr_mean, times = 100)
# resample1



null_dist <- gss %>%
    specify(response = hours) %>%
    hypothesize(null = "point", mu = 40) %>%
    generate(reps = 5000, type = "bootstrap") %>%
    calculate(stat = "mean")

point_estimate <- gss %>%
    specify(response = hours) %>%
    calculate(stat = "mean")

null_dist %>%
    get_p_value(obs_stat = point_estimate, direction = "two_sided")




null_dist2 <- lr_cpdb2 %>%
    unite(col = "cell_pair", source, target, sep = "~") %>%
    group_by(cell_pair)

null_dist2 %>%
    group_split(.keep = TRUE) %>%
    map(function(x)
        x %>% specify(response = lr_mean) %>%
            hypothesize(null = "point", mu = mean(lr_cpdb2[["lr_mean"]])) %>%
            generate(reps = 100, type = "bootstrap") %>%
            calculate(stat = "mean"))  %>%
    setNames(unique(null_dist2$cell_pair))




point_estimate <- gss %>%
    specify(response = hours) %>%
    calculate(stat = "mean")




point_estimate2 <- lr_cpdb2 %>%
    rowwise() %>%
    specify(response = lr_mean) %>%
    calculate(stat = "mean")


xd <- lr_cpdb2 %>%
    mutate(lr_pval =
               null_dist2 %>%
               get_p_value(obs_stat = lr_mean, direction = "two_sided") %>%
               pull(p_value))

null_dist2 %>%
    get_p_value(obs_stat = point_estimate2, direction = "two_sided")


