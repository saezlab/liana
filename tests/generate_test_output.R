# input
complex_policy='mean0'
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                                   "input", "testdata.rds"))

# convert and save the Seurat object
seurat_conv <- liana_prep(seurat_object)
# saveRDS(seurat_conv, file.path(liana_path , "testdata",
#                                "input", "seurat_conv.rds"))

# # Generate SingleCellExperiment Object
# set.seed(123)
# counts <- matrix(replicate(90,rpois(100, lambda = 10)), ncol=90, nrow=1000)
# pretend.cell.labels <- colnames(seurat_object)
# pretend.gene.lengths <- sample(10000, nrow(counts))
#
# sce <- SingleCellExperiment(list(counts=counts),
#                             colData=data.frame(label=pretend.cell.labels),
#                             rowData=data.frame(length=pretend.gene.lengths),
#                             metadata=list(study="GSE111111")
# )
# rownames(sce) <- rownames(seurat_object)[1:nrow(counts)]
# colnames(sce) <- pretend.cell.labels
# sce <- scater::logNormCounts(sce)
# colLabels(sce) <- factor(sample(letters[1:3], ncol(counts), replace=TRUE))
# saveRDS(sce, file.path(liana_path , "testdata",
#                        "input", "testsce.rds"))

# liana PREP output ----
# test sce as input
sce <- readRDS(file.path(liana_path , "testdata",
                         "input", "testsce.rds"))

sce_conv <- liana_prep(sce)
saveRDS(sce_conv, file.path(liana_path , "testdata",
                            "input", "sce_conv.rds"))

wrap_sce <- liana_wrap(sce, resource = "OmniPath", method = c("sca", "natmi"))
saveRDS(wrap_sce, file.path(liana_path , "testdata",
                            "input", "wrap_sce.rds"))


# liana Pipe Output ----
pipe_out <- liana_pipe(seurat_conv,
                       op_resource = select_resource("OmniPath")[[1]] %>%
                           decomplexify())
saveRDS(pipe_out, file.path(liana_path, "testdata",
                            "output", "liana_pipe.RDS"))


# Scores Output ----
conn_out <- get_connectome(pipe_out, expr_prop=0, complex_policy=complex_policy)
saveRDS(conn_out, file.path(liana_path, "testdata",
                            "output", "conn_score.RDS"))

logfc_out <- get_logfc(pipe_out, complex_policy=complex_policy)
saveRDS(logfc_out, file.path(liana_path, "testdata",
                            "output", "logfc_score.RDS"))

natmi_out <- get_natmi(pipe_out, complex_policy=complex_policy)
saveRDS(natmi_out, file.path(liana_path, "testdata",
                            "output", "natmi_score.RDS"))

sca_out <- get_sca(pipe_out, complex_policy=complex_policy)
saveRDS(sca_out, file.path(liana_path, "testdata",
                           "output", "sca_score.RDS"))

# liana permutations and cpdb output
cpdb_out <- liana_wrap(seurat_object,
                       method = c('cellphonedb'),
                       resource = c('CellPhoneDB'),
                       permutation.params = list(nperms=20))
saveRDS(cpdb_out, file.path(liana_path, "testdata",
                            "output", "liana_cpdb.RDS"))

# cytotalk
cytotalk_out <- liana_wrap(seurat_object,
                           method = c('cytotalk'),
                           resource = c('OmniPath'))
saveRDS(cytotalk_out, file.path(liana_path, "testdata",
                                "output", "liana_cytotalk.RDS"))

# Recomplexify Output ----
lr_cmplx <- liana_pipe(seurat_conv,
                       op_resource =
                           select_resource("CellPhoneDB")[[1]] %>%
                           decomplexify())

recomplex <- recomplexify(lr_cmplx,
                          .score_specs()[["sca"]]@columns,
                          complex_policy =complex_policy)
saveRDS(recomplex, file.path(liana_path, "testdata",
                             "output", "recomplex.RDS"))

# liana Wrapper Output ----
wrap_out <- liana_wrap(seurat_object,
                       method = c('logfc','natmi', 'connectome'),
                       resource = c('OmniPath'))
saveRDS(wrap_out, file.path(liana_path, "testdata",
                             "output", "liana_res.RDS"))

# wrap_default
wrap_def_out <- liana_wrap(seurat_object,
                           method = c('sca','squidpy', "call_sca"),
                           resource = "Default")
saveRDS(wrap_def_out, file.path(liana_path, "testdata",
                                "output", "liana_def_res.RDS"))


# LIANA Defaults
def_arg <- liana_defaults(expr_prop=0,
                          squidpy.params=list(threshold = 0.1),
                          cellchat.params=list(nboot=1000))
saveRDS(def_arg, file.path(liana_path, "testdata",
                           "output", "liana_def_args.RDS"))

# liana dotplot ----
liana_dotplot_out <- liana_dotplot(cpdb_out %>%
                                       mutate(pvalue = -log10(pvalue+0.00000001)), # invert pval
                                   source_groups = "B",
                                   target_groups = c("NK", "CD8 T"),
                                   magnitude = "lr.mean",
                                   specificity = "pvalue",
                                   show_complex = TRUE)
saveRDS(liana_dotplot_out,
        file.path(liana_path, "testdata",
                  "output", "liana_dotplot_out.RDS"))


### EXTERNAL ----
# Test Connectome ----
conn_res <- call_connectome(
    sce = seurat_object,
    op_resource = select_resource("OmniPath")[[1]], # Default = No sig hits
    .spatial = FALSE,
    min.cells.per.ident = 1,
    p.values = TRUE,
    calculate.DOR = FALSE,
    assay = 'RNA',
    .format = TRUE
)
saveRDS(conn_res, file.path(liana_path, "testdata",
                            "output", "conn_res.RDS"))

# Test Squidpy ----
squidpy_res <- call_squidpy(sce = seurat_object,
                            op_resource = select_resource("OmniPath"),
                            cluster_key="seurat_annotations",
                            n_perms=100,
                            threshold=0.01,
                            seed=as.integer(1004))
saveRDS(squidpy_res, file.path(liana_path, "testdata",
                            "output", "squidpy_res.RDS"))


# Test NATMI ----
natmi_res <- call_natmi(op_resource = select_resource("OmniPath")[[1]],
                        sce = seurat_object,
                        expr_file = "test_em.csv",
                        meta_file = "test_metadata.csv",
                        output_dir = "NATMI_test",
                        assay = "RNA",
                        num_cor = 4,
                        .format = TRUE,
                        .overwrite_data = TRUE,
                        .seed = 1004,
                        .natmi_path = NULL)
saveRDS(natmi_res, file.path(liana_path, "testdata",
                             "output", "natmi_res.RDS"))

# Test iTALK ----
italk_res <- call_italk(op_resource = NULL,
                        sce = seurat_object,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)
saveRDS(italk_res, file.path(liana_path, "testdata",
                             "output", "italk_res.RDS"))


# Test SCA ----
sca_res <- call_sca(op_resource = NULL,
                    sce = seurat_object,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5))
saveRDS(sca_res, file.path(liana_path, "testdata",
                           "output", "sca_res.RDS"))


# Test CellChat ----
cellchat_res <- call_cellchat(
    op_resource = NULL,
    sce = seurat_object,
    nboot = 2,
    exclude_anns = NULL,
    thresh = 1,
    assay = "RNA",
    .normalize = FALSE,
    .do_parallel = FALSE,
    .raw_use = TRUE)

saveRDS(cellchat_res,
        file.path(liana_path, "testdata",
                  "output", "cc_res.RDS"))

# liana aggregate output ----
# Simplest scenario
liana_aggr <- readRDS(file.path(liana_path, "testdata",
                                "output", "liana_res.RDS")) %>%
    liana_aggregate()
saveRDS(liana_aggr, file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))

# External Methods + Housekeep scores
liana_res <- readRDS(file.path(liana_path, "testdata",
                               "output", "liana_res.RDS"))

liana_res$call_natmi <- readRDS(file.path(liana_path, "testdata",
                                          "output", "natmi_res.RDS"))

liana_res$cellchat <- readRDS(file.path(liana_path, "testdata",
                                        "output", "cc_res.RDS"))
saveRDS(liana_res, file.path(liana_path, "testdata",
                             "output", "liana_res_plus.RDS"))

liana_agg_house <- liana_res %>%
    liana_aggregate(.score_mode = .score_housekeep)
saveRDS(liana_agg_house, file.path(liana_path, "testdata",
                              "output", "liana_house_aggr.RDS"))


