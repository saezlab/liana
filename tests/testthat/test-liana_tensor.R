# Input
liana_path <- system.file(package = "liana")
context_df_dict <- readRDS(file.path(liana_path, "testdata",
                                     "output", "context_df_test.RDS"))

test_that("Test tensor wrapper", {
    skip_on_ci()
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "decomp.RDS"))
    res1 <- liana_tensor_c2c(context_df_dict = context_df_dict,
                             score_col = 'score',
                             rank = 5,
                             how='outer',
                             inplace = FALSE
    )

    sum_res1 <- round(colSums(res1$contexts[, sapply(res1$contexts, is.numeric)]), 5)
    sum_exp1 <- round(colSums(exp1$contexts[, sapply(exp1$contexts, is.numeric)]), 5)
    expect_true(identical(sum_res1, sum_exp1))

})

test_that("Test decompose_tensor", {
    tensor <- liana_tensor_c2c(context_df_dict = context_df_dict,
                               score_col = 'score',
                               rank = 5,
                               how='outer',
                               build_only = TRUE,
                               inplace = FALSE
    )
    res1 <- decompose_tensor(tensor =  tensor,
                             tf_optimization='regular',
                             upper_rank = 5)
    expect_equal(dim(res1$contexts), c(2, 2))
})

