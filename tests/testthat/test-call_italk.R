# # Input
# liana_path <- system.file(package = "liana")
# seurat_object <-
#     readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
#
# # Test
# test_that("Test iTALK", {
#     if(!exists("test_external")){
#         skip("Not testing externals")
#     } else{
#         if(!test_external) skip("Not testing externals")
#     }
#
#     exp1 <- readRDS(file.path(liana_path, "testdata",
#                               "output", "italk_res.RDS"))
#     res1 <- suppressWarnings(call_italk(sce = seurat_object,
#                                         op_resource = NULL,
#                                         assay = 'RNA',
#                                         .format = TRUE,
#                                         .DE = TRUE))
#     expect_equal(exp1, res1)
# })
