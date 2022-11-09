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

    # fix random bs exception
    exp1$interactions$lr <- as.character(exp1$interactions$lr)
    res1$interactions$lr <- as.character(res1$interactions$lr)

    expect_true(all(map2_lgl(res1, exp1, ~all_equal(.x, .y))))
})
