# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# liana Pipe Output ----
pipe_out <- liana_pipe(seurat_object,
                       op_resource = select_resource("OmniPath")[[1]])
saveRDS(pipe_out, file.path(liana_path, "testdata",
                            "output", "liana_pipe.RDS"))


# Scores Output ----
conn_score <- get_connectome(pipe_out)
saveRDS(conn_score, file.path(liana_path, "testdata",
                            "output", "conn_score.RDS"))

logfc_score <- get_logfc(pipe_out)
saveRDS(logfc_score, file.path(liana_path, "testdata",
                            "output", "logfc_score.RDS"))

natmi_score <- get_natmi(pipe_out)
saveRDS(natmi_score, file.path(liana_path, "testdata",
                            "output", "natmi_score.RDS"))

sca_score <- get_sca(pipe_out)
saveRDS(sca_score, file.path(liana_path, "testdata",
                              "output", "sca_score.RDS"))



# Recomplexify Output ----
recomplex <- recomplexify(pipe_out,
                          .score_specs()[["sca"]]@columns,
                          protein = 'subunit',
                          complex_policy ='min0')
saveRDS(recomplex, file.path(liana_path, "testdata",
                             "output", "recomplex.RDS"))

# liana Wrapper Output ----
wrap_out <- liana_wrap(seurat_object,
                        method = c('sca','squidpy'),
                        resource = c('OmniPath'))
saveRDS(wrap_out, file.path(liana_path, "testdata",
                             "output", "liana_res.RDS"))

# liana aggregate output ----
liana_aggr <- readRDS(file.path(liana_path, "testdata",
                          "output", "liana_res.RDS")) %>%
    liana_aggregate()
saveRDS(liana_aggr, file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))

# Test Connectome ----
conn_res <- call_connectome(
    seurat_object = seurat_object,
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
squidpy_res <- call_squidpy(seurat_object = seurat_object,
                            op_resource = select_resource("OmniPath"),
                            cluster_key="seurat_annotations",
                            n_perms=100,
                            threshold=0.01,
                            seed=as.integer(1004))
saveRDS(squidpy_res, file.path(liana_path, "testdata",
                            "output", "squidpy_res.RDS"))


# Test NATMI ----
natmi_res <- call_natmi(op_resource = select_resource("OmniPath"),
                        seurat_object = seurat_object,
                        expr_file = "test_em.csv",
                        meta_file = "test_metadata.csv",
                        output_dir = "NATMI_test",
                        assay = "RNA",
                        num_cor = 4,
                        .format = TRUE,
                        .write_data = TRUE,
                        .seed = 1004,
                        .natmi_path = NULL)
saveRDS(natmi_res, file.path(liana_path, "testdata",
                             "output", "natmi_res.RDS"))

# Test iTALK ----
italk_res <- call_italk(op_resource = NULL,
                        seurat_object = seurat_object,
                        assay = 'RNA',
                        .format = TRUE,
                        .DE = TRUE)
saveRDS(italk_res, file.path(liana_path, "testdata",
                             "output", "italk_res.RDS"))


# Test SCA ----
sca_res <- call_sca(op_resource = NULL,
                    seurat_object = seurat_object,
                    assay = 'RNA',
                    .format = TRUE,
                    s.score = 0,
                    logFC = log2(1.5))
saveRDS(sca_res, file.path(liana_path, "testdata",
                           "output", "sca_res.RDS"))

