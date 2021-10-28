# Input
liana_path <- system.file(package = "liana")

# Test without decomplexify
test_that("Test liana aggregate", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_aggr.RDS"))

    res1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS")) %>%
        liana_aggregate(.decomplexify = FALSE)

    expect_equal(exp2, res2)
})


# Test with decomplexify
test_that("Test liana aggregate with CellChat complexes", {
    res2 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS"))
    # append cellchat
    res2$cellchat <- readRDS(file.path(liana_path, "testdata",
                                       "output", "cc_res.RDS"))

    res2_no <- res2 %>% liana_aggregate(.decomplexify = FALSE)

    res2 %<>% liana_aggregate()

    expect_equal(exp1, res1)
})
