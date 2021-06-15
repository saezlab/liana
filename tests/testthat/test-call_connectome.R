# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test Connectome", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "conn_res.RDS"))
    res1 <- call_connectome(
        seurat_object = seurat_object,
        op_resource = select_resource("OmniPath")[[1]], # Default = No sig hits
        .spatial = FALSE,
        min.cells.per.ident = 1,
        p.values = TRUE,
        calculate.DOR = FALSE,
        assay = 'RNA',
        .format = TRUE
    )

    expect_equal(exp1, res1)
})
