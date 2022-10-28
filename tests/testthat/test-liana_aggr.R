# Input
liana_path <- system.file(package = "liana")

# Test without decomplexify ----
test_that("Test liana aggregate", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))

    res1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS")) %>%
        liana_aggregate()

    expect_equal(exp1, res1)
})


# aggregate with housekeeping ----
test_that("Test liana aggregate (housekeep)", {
    exp2 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_house_aggr.RDS"))

    # warning should be thrown by liana-sca not having a 'housekeeping score'
    testthat::expect_warning(
        res2 <- readRDS(file.path(liana_path, "testdata",
                                  "output", "liana_res_plus.RDS")) %>%
            liana_aggregate(.score_mode = .score_housekeep)
    )


    expect_equal(exp2, res2)
})



# Rank method test ----
test_that("Test Rank Method", {

    liana_res <- readRDS(file.path(liana_path, "testdata",
                                   "output", "liana_res.RDS"))

    spec <- readRDS(file.path(liana_path, "testdata",
                      "output", "rank_spec_example.RDS"))
    spec_test <- liana_res$natmi %>%
        rank_method("natmi", "specificity")
    mag <- readRDS(file.path(liana_path, "testdata",
                              "output", "rank_mag_example.RDS"))
    mag_test <- liana_res$natmi %>%
        rank_method("natmi", "magnitude")


    expect_equal(spec, spec_test)
    expect_equal(mag, mag_test)
    expect_error(rank_method(liana_res$natmi,"natmi", "xd"),
                 "Passed `mode` not found!")
    expect_error(rank_method(liana_res$logfc,"logfc", "magnitude"),
                 "Score-method combination not found!")
})

