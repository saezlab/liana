# Input
liana_path <- system.file(package = "liana")

# Test
test_that("Test liana aggregate", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))
    res1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS")) %>%
        liana_aggregate()

    expect_equal(exp1, res1)
})
