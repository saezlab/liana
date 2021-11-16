# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test with OmniPath ----
test_that("Test liana wrapper", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS"))
    res1 <- liana_wrap(seurat_object,
                       method = c('logfc','natmi', 'connectome'),
                       resource = c('OmniPath'))

    expect_equal(res1, exp1)
})

# Test with Default ----
test_that("Test liana wrapper", {
    exp2 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_def_res.RDS"))
    res2 <- liana_wrap(seurat_object,
                       method = c('sca','squidpy', "call_sca"),
                       resource = "Default")

    expect_equal(exp2, res2)
})


# Test Default Params
test_that("Test liana wrapper", {
    exp3 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_def_args.RDS"))
    res3 <- liana_defaults(expr_prop=0,
                           squidpy.params=list(threshold = 0.1),
                           cellchat.params=list(nboot=1000))

    expect_equal(exp3, res3)
})

# Test Expression_prop filtering
test_that("Test expr_prop filtering", {
    liana_pipe_res <- readRDS(file.path(liana_path, "testdata",
                                        "output", "liana_pipe.RDS"))

    expect_equal(10548, nrow(.filt_liana_pipe(liana_pipe_res, "connectome")))
    expect_equal(998, nrow(.filt_liana_pipe(liana_pipe_res, "cellphonedb")))

})


