# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test liana wrapper", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS"))
    res1 <- liana_wrap(seurat_object,
                       method = c('sca','squidpy'),
                       resource = c('OmniPath'))

    expect_equal(exp1, res1)
})
