# Input
liana_path <- system.file(package = "liana")

# Test without decomplexify
test_that("Test liana aggregate", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))

    res1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS")) %>%
        liana_aggregate()

    expect_equal(exp1, res1)
})


# aggregate with housekeeping
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
