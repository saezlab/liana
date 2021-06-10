# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

omni_resources <- readRDS(file.path(liana_path , "omni_resources.rds"))
op_resource <- omni_resources["OmniPath"]

# Test
test_that("Test NATMI", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "natmi_res.RDS"))
    res1 <- call_natmi(op_resource = op_resource,
                       seurat_object = seurat_object,
                       expr_file = "em.csv",
                       meta_file = "metadata.csv",
                       output_dir = "NATMI_test",
                       assay = "RNA",
                       num_cor = 4,
                       .format = TRUE,
                       .write_data = FALSE,
                       .seed = 1004,
                       .natmi_path = NULL)

    expect_equal(exp1, res1)
})
