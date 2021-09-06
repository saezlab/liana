# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds")) %>%
    Seurat::NormalizeData()
op_resource <- select_resource("CellChatDB")[[1]]

# Resource Format
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
entity_genes = union(transmitters$gene,
                     receivers$gene)

score_col = "pvalue"
parallelize = F
workers = 4

## CPDB - this part is called in liana_wrap
lr_res <- liana_pipe(seurat_object = seurat_object,
                     op_resource = op_resource %>% decomplexify(),
                     expr_prop = 0.1,
                     trim = 0,
                     assay.type = "logcounts")

sce <- seurat_to_sce(seurat_object = seurat_object,
                     entity_genes,
                     assay="RNA")


# reproduce liana_scores
score_object <- .score_specs()[["cellphonedb"]]
complex_policy = "min0"

lr_res %<>%
    select(ligand, receptor,
           ends_with("complex"),
           source, target,
           !!score_object@columns)

lr_res %<>%
    recomplexify(
        lr_res = .,
        columns = score_object@columns,
        complex_policy = complex_policy)


perm_means <- get_permutations(lr_res,
                               sce,
                               nperms=100,
                               seed=1234,
                               trim=0.1,
                               parallelize = FALSE,
                               workers=4)


dotdotdot <- list(parallelize = FALSE,
                  workers = 4,
                  perm_means = perm_means)

args <-
    append(
        list(lr_res = lr_res,
             score_col = score_object@method_score),
        dotdotdot
    )

liana_cpdb <- exec(
    cellphonedb_score,
    !!!args) %>%
    ungroup()



# Check perms
parallelize = FALSE
workers = 4

og_res <- lr_res %>%
    rowwise() %>%
    select(ligand, receptor, source, target)

progress_bar <- dplyr::progress_estimated(length(perm_means) * 2)
perm_joined <- perm_means %>%
    map_custom(function(pmean){
        og_res %>%
            distinct() %>%
            liana:::join_means(means = pmean,
                               source_target = "source",
                               entity = "ligand",
                               type = "trunc",
                               pb = progress_bar) %>%
            liana:::join_means(means = pmean,
                               source_target = "target",
                               entity = "receptor",
                               type = "trunc",
                               pb = progress_bar) %>%
            replace(is.na(.), 0) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(lr_mean = mean(c(ligand.trunc, receptor.trunc)))
    }, parallelize = parallelize, workers = workers) %>%
    bind_rows()

#
quantiles <- lr_res %>%
    group_by(source, target, ligand.complex, receptor.complex) %>%
    mutate(og_mean = mean(c(ligand.trunc, receptor.trunc))) %>%
    select(source, target,
           ligand.complex, ligand,
           receptor, receptor.complex,
           og_mean) %>%
    left_join(perm_joined, by = c("source", "target", "ligand", "receptor")) %>%
    group_by(ligand.complex, receptor.complex, source, target) %>%
    group_split() %>%
    map_dbl(function(interaction){
        null_dist <- ecdf(interaction$lr_mean)
        og_mean <- interaction %>% pull("og_mean") %>% unique()
        null_dist(og_mean) %>% pluck(1)
    })


# quantiles <- perm_joined %>%
#     group_by(ligand, receptor, source, target) %>%
#     group_split() %>%
#     map_dbl(function(interaction){
#         null_dist <- ecdf(interaction$lr_mean)
#         og_mean <- interaction %>% pull("og_mean") %>% unique()
#         null_dist(og_mean) %>% pluck(1)
#     })
#



pvals_df <- lr_res %>%
    group_by(ligand.complex, receptor.complex, source, target) %>%
    group_keys() %>%
    mutate( {{ score_col }} :=  1 - quantiles)


liana_cpdb4 <- lr_res %>%
    rowwise() %>%
    mutate(lr.mean = mean(c(ligand.trunc, receptor.trunc))) %>%
    left_join(pvals_df,
              by = c("ligand.complex", "receptor.complex",
                     "source", "target")) %>%
    mutate({{ score_col }} := # replace pval of non-expressed rec and ligs
               ifelse(ligand.trunc == 0 || receptor.trunc == 0,
                      1,
                      .data[[score_col]]))



perm_joined %>%
    filter(ligand == "TGFB1" && receptor == "TGFBR2" &&
               source == "CD8 T" && target == "CD8 T") %>%
    filter(lr_mean >= 0.0774)







liana_cpdb4 %>%
    filter(ligand == "TGFB1" && receptor == "TGFBR1" &&
               source == "NK" && target == "NK")

perm_joined %>%
    filter(ligand == "TGFB1" && receptor == "TGFBR1" &&
               source == "NK" && target == "NK") %>%
    filter(lr_mean >= 0.455)



### via liana_call
liana_cpdb2 <- liana_call("cellphonedb",
                          lr_res = lr_res,
                          seurat_object = seurat_object,
                          perm_means = perm_means,
                          parallelize = FALSE,
                          workers = 4)


### via liana_wrap
liana_cpdb3 <- liana_wrap(seurat_object = seurat_object,
                          method = "cellphonedb",
                          resource = "CellPhoneDB",
                          permutation.params = list(nperms=100),
                          expr_prop = 0.1
                          )


#
op_resources <- select_resource("OmniPath")[[1]] %>%
                        select(
                            uniprot_source = source,
                            unprot_target = target,
                            source = source_genesymbol,
                            target = target_genesymbol,
                            category_intercell_source,
                            category_intercell_target
                        )
write.csv(op_resources, "../cpdb/input/op_resource.csv", row.names = FALSE)


# check
lr_res

liana_cpdb3 <- liana_wrap(seurat_object = seurat_object,
                          method = c("cellphonedb", "squidpy"),
                          resource = "CellPhoneDB",
                          permutation.params = list(nperms=1000),
                          expr_prop = 0.1,
                          trim = 0)

xd <- liana_cpdb3 %>%
    liana_aggregate()


liana_cpdb3$cellphonedb %>% filter(pvalue <= 0.05)
liana_cpdb3$squidpy %>% filter(pvalue <= 0.05)
og_cpdb %>% filter(pvalue <= 0.05)


mean_prop <-
    scuttle::summarizeAssayByGroup(sce,
                                   ids = colLabels(sce),
                                   assay.type = "logcounts",
                                   statistics = c("mean", "prop.detected"))
means <- mean_prop@assays@data$mean %>%
    as_tibble(rownames = "gene")
props <- mean_prop@assays@data$prop.detected %>%
    as_tibble(rownames = "gene")



#
liana_cpdb3 <- liana_wrap(seurat_object = seurat_object,
                          method = c("cellphonedb", "squidpy"),
                          resource = "CellPhoneDB",
                          permutation.params = list(nperms=1000),
                          expr_prop = 0,
                          trim = 0,
                          assay.type = "counts")
liana_cpdb3$cellphonedb %>% filter(pvalue <= 0.05)
liana_cpdb3$squidpy %>% filter(pvalue <= 0.05)

xd <- liana_cpdb3 %>% liana_aggregate()
