# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

omni_resources <- readRDS(file.path(liana_path , "omni_resources.rds"))
op_resource <- omni_resources["OmniPath"]

# Test
test_that("Test Squidpy", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "squidpy_res.RDS"))
    res1 <- call_squidpyR(seurat_object = seurat_object,
                          op_resource = op_resource,
                          cluster_key="seurat_annotations",
                          n_perms=100,
                          threshold=0.01,
                          seed=as.integer(1004))

    expect_equal(exp1, res1)
})
