# Input----
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
pipe_out <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_pipe.RDS"))

# Test Liana Pipe----
test_that("Test liana pipe", {
    res1 <- liana_pipe(sce = liana_prep(sce = seurat_object),
                       op_resource = select_resource("OmniPath")[[1]] %>%
                           decomplexify(),
                       base=2)

    expect_equal(res1, pipe_out)
})


# Test De-/Re- Complexify----
test_that("Test liana pipe", {
    lr_cmplx <- liana_pipe(liana_prep(sce = seurat_object),
                           op_resource = select_resource("CellPhoneDB")[[1]] %>%
                               decomplexify(),
                           base=2)

    recomplex_exp <- readRDS(file.path(liana_path, "testdata",
                                       "output", "recomplex.RDS"))
    recomplex <- recomplexify(lr_cmplx,
                              .score_specs()[["sca"]]@columns,
                              complex_policy ='mean0')

    expect_equal(recomplex, recomplex_exp)
})


# Test Get Scores
test_that("Test LIANA Scores", {
    complex_policy='mean0'

    conn_score <- get_connectome(pipe_out,
                                 complex_policy=complex_policy,
                                 expr_prop=liana_defaults()[["expr_prop"]])
    conn_exp <- readRDS(file.path(liana_path, "testdata",
                                  "output", "conn_score.RDS"))

    logfc_score <- get_logfc(pipe_out,
                             complex_policy=complex_policy,
                             expr_prop=liana_defaults()[["expr_prop"]])
    logfc_exp <- readRDS(file.path(liana_path, "testdata",
                                   "output", "logfc_score.RDS"))

    natmi_score <- get_natmi(pipe_out,
                             complex_policy=complex_policy,
                             expr_prop=liana_defaults()[["expr_prop"]])
    natmi_exp <- readRDS(file.path(liana_path, "testdata",
                                   "output", "natmi_score.RDS"))

    sca_score <- get_sca(pipe_out,
                         complex_policy=complex_policy,
                         expr_prop=liana_defaults()[["expr_prop"]])
    sca_exp <- readRDS(file.path(liana_path, "testdata",
                                 "output", "sca_score.RDS"))

    expect_equal(conn_score, conn_exp)
    expect_equal(logfc_score, logfc_exp)
    expect_equal(natmi_score, natmi_exp)
    expect_equal(sca_score, sca_exp)
})

