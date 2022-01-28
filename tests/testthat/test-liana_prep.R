# Input
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                                   "input", "testdata.rds"))
sce  <- readRDS(file.path(liana_path , "testdata",
                          "input", "testsce.rds"))

# Test Conversion and filter
test_that("Test liana with SingleCellExperiment", {
    exp1 <- readRDS(file.path(liana_path , "testdata",
                              "input", "sce_conv.rds"))
    res1 <- liana_prep(sce)

    exp2 <- readRDS(file.path(liana_path , "testdata",
                              "input", "wrap_sce.rds"))
    res2 <- liana_wrap(sce, resource = "OmniPath", method = c("sca", "natmi"))

    suppressMessages(all_equal(res1@assays@data$counts, exp1@assays@data$counts))
    all_equal(res2[[1]], exp2[[1]])
    all_equal(res2[[2]], exp2[[2]])
})

test_that("Test Seurat to SCE conversion", {
    res3 <- liana_prep(seurat_object)
    exp3 <- readRDS(file.path(liana_path , "testdata",
                                   "input", "seurat_conv.rds"))

    all_equal(res3@assays@data$logcounts,
              exp3@assays@data$logcounts)
    expect_warning(liana_prep(seurat_object))
})
