# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
pipe_out <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_pipe.RDS"))

# Test Liana Pipe
test_that("Test liana pipe", {
    res1 <- liana_pipe(seurat_object = seurat_object,
                       op_resource = select_resource("OmniPath")[[1]] %>%
                           decomplexify())

    expect_equal(pipe_out, res1)
})

# Test De-/Re- Complexify
test_that("Test liana pipe", {
    lr_cmplx <- liana_pipe(seurat_object,
                           op_resource = select_resource("CellPhoneDB")[[1]] %>%
                               decomplexify())

    recomplex_exp <- readRDS(file.path(liana_path, "testdata",
                                       "output", "recomplex.RDS"))
    recomplex <- recomplexify(lr_cmplx,
                              .score_specs()[["sca"]]@columns,
                              complex_policy ='min0')

    expect_equal(recomplex_exp, recomplex)
})



# Test Get Scores
test_that("Test LIANA Scores", {
    conn_score <- get_connectome(pipe_out, expr_prop=0)
    conn_exp <- readRDS(file.path(liana_path, "testdata",
                                  "output", "conn_score.RDS"))

    logfc_score <- get_logfc(pipe_out)
    logfc_exp <- readRDS(file.path(liana_path, "testdata",
                                   "output", "logfc_score.RDS"))

    natmi_score <- get_natmi(pipe_out)
    natmi_exp <- readRDS(file.path(liana_path, "testdata",
                                   "output", "natmi_score.RDS"))

    sca_score <- get_sca(pipe_out)
    sca_exp <- readRDS(file.path(liana_path, "testdata",
                                 "output", "sca_score.RDS"))

    expect_equal(conn_score, conn_exp)
    expect_equal(logfc_score, logfc_exp)
    expect_equal(natmi_score, natmi_exp)
    expect_equal(sca_score, sca_exp)
})

