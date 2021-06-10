# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("multiplication works", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "cc_res.RDS"))
    res1 <- call_cellchat(
        op_resource = NULL,
        seurat_object = seurat_object,
        nboot = 10,
        exclude_anns = NULL,
        thresh = 1,
        assay = "RNA",
        .normalize = TRUE,
        .do_parallel = FALSE,
        .raw_use = TRUE
    )

    expect_equal(exp1, res1)
})
