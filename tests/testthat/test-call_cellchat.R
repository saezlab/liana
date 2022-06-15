# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))


test_that("Test CellChat with Default", {
    if(!exists("test_external")){
        skip("Not testing externals")
    } else{
        if(!test_external) skip("Not testing externals")
    }
    skip_on_ci()

    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "cc_res.RDS"))
    res1 <- call_cellchat(
        op_resource = NULL,
        sce = seurat_object,
        nboot = 2,
        exclude_anns = NULL,
        thresh = 1.1, # cellchat filtering changed as of 1.1.3, I assume to <
        assay = "RNA",
        .normalize = FALSE,
        .do_parallel = FALSE,
        .raw_use = TRUE
    )

    expect_equal(res1, exp1)

})
