# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test
test_that("Test NATMI", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "natmi_res.RDS"))
    res1 <- call_natmi(op_resource = select_resource("OmniPath")[[1]],
                       seurat_object = seurat_object,
                       expr_file = "test_em.csv",
                       meta_file = "test_metadata.csv",
                       output_dir = "NATMI_test",
                       assay = "RNA",
                       num_cor = 4,
                       .format = TRUE,
                       .seed = 1004,
                       .natmi_path = NULL,
                       .overwrite_data = FALSE)

    expect_equal(exp1, res1)
})


