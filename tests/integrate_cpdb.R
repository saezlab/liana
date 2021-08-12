# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds")) %>%
    Seurat::NormalizeData()
op_resource <- select_resource("CellPhoneDB")[[1]]
# Resource Format
transmitters <- op_resource$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
receivers <- op_resource$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
entity_genes = union(transmitters$gene,
                     receivers$gene)


## CPDB - this part is called in liana_wrap
lr_res <- liana_pipe(seurat_object = seurat_object,
                     op_resource = op_resource %>% decomplexify(),
                     expr_prop = 0.2,
                     trim = 0.1,
                     assay.type = "logcounts")

sce <- seurat_to_sce(seurat_object = seurat_object,
                     entity_genes,
                     assay="RNA")


# Run alg externally: - this does not deal with complexes
perm_means_ext <- get_permutations(lr_res,
                                   sce,
                                   nperms=10,
                                   seed=1234,
                                   trim=0.1,
                                   parallelize = FALSE,
                                   workers=4)

liana_cpdb_ext <- cellphonedb_score(lr_res = lr_res,
                                    perm_means = perm_means_ext,
                                    parallelize = FALSE,
                                    workers = 4,
                                    score_col = "pvalue")




# reproduce liana_scores
#
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


### via liana_call
xx <- liana_call("cellphonedb",
                 lr_res = lr_res,
                 seurat_object = seurat_object,
                 perm_means = perm_means,
                 parallelize = FALSE,
                 workers = 4)


### via liana_wrap
liana_cpdb <- liana_wrap(seurat_object = seurat_object,
                         method = "cellphonedb")


