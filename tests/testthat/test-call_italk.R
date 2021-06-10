# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test iTALK", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "italk_res.RDS"))
    res1 <- call_italk(op_resource = NULL,
                       seurat_object = seurat_object,
                       assay = 'RNA',
                       .format = TRUE,
                       .DE = TRUE)

    expect_equal(exp1, res1)
})
