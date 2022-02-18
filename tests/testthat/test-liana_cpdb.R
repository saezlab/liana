# Input
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                                   "input", "testdata.rds"))

# Test
test_that("Test liana cpdb", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_cpdb.RDS"))
    res1 <- liana_wrap(seurat_object,
                       method = c('cellphonedb'),
                       resource = c('Consensus'),
                       permutation.params = list(nperms=20))

    expect_equal(exp1, res1)
})
