# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test SingleCellSignalR", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "sca_res.RDS"))
    res1 <- call_sca(op_resource = NULL,
                     seurat_object = seurat_object,
                     assay = 'RNA',
                     .format = TRUE,
                     s.score = 0,
                     logFC = log2(1.5)
    )

    expect_equal(exp1, res1)
})

