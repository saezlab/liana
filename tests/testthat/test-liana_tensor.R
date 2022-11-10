# Input
liana_path <- system.file(package = "liana")
context_df_dict <- readRDS(file.path(liana_path, "testdata",
                                     "output", "context_df_test.RDS"))

test_that("Test tensor wrapper", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "decomp.RDS"))
    res1 <- liana_tensor_c2c(context_df_dict = context_df_dict,
                             score_col = 'score',
                             rank = 5,
                             how='outer',
                             inplace = FALSE
                             )

    # fix random bs exceptions
    exp1 <- map(exp1, ~.x %>% mutate(across(where(is.factor), ~as.character(.x))))
    res1 <- map(res1, ~.x %>% mutate(across(where(is.factor), ~as.character(.x))))

    expect_true(all_equal(res1$contexts, exp1$contexts))
    expect_true(all_equal(res1$interactions, exp1$interactions))
    expect_true(all_equal(res1$senders, exp1$senders))
    expect_true(all_equal(res1$receivers, exp1$receivers))
})
