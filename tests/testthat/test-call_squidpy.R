# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test Squidpy", {
    if(!exists("test_external")){
        skip("Not testing externals")
    } else{
        if(!test_external) skip("Not testing externals")
    }
    skip_on_ci()

    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "squidpy_res.RDS"))
    res1 <- call_squidpy(sce = seurat_object,
                         op_resource = select_resource("OmniPath")[[1]],
                         cluster_key=NULL,
                         n_perms=100,
                         threshold=0.01,
                         seed=as.integer(1004))

    expect_equal(exp1, res1)
})
