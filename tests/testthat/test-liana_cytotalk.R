liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))

# Test Cytotalk Wrap
test_that("Test Cytotalk Wrap", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_cytotalk.RDS"))

    res1 <- liana_wrap(seurat_object,
                       method = c('cytotalk'),
                       resource = c('OmniPath'))

    expect_equal(exp1, res1)
})
