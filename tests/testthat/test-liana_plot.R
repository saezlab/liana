# Input
liana_path <- system.file(package = "liana")
liana_cpdb <- readRDS(file.path(liana_path, "testdata",
                                "output", "liana_cpdb.RDS"))
liana_aggr <- readRDS(file.path(liana_path, "testdata",
                                "output", "liana_aggr.RDS")) %>%
    filter(aggregate_rank <= 0.05)

# Test DotPlot
test_that("Test liana dotplot", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_dotplot_out.RDS"))
    res1 <- liana_dotplot(liana_cpdb %>% # invert
                              mutate(pvalue = -log10(pvalue+1e10)),
                          source_groups = "B",
                          target_groups = c("NK", "CD8 T"),
                          magnitude = "lr.mean",
                          specificity = "pvalue",
                          show_complex = TRUE)

    expect_identical(res1$data, exp1$data)
})

# Test Heatmaps
test_that("Test liana heatmaps", {
    exp_freq <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_freq_out.RDS"))

    res_freq <- heat_freq(liana_cpdb %>% filter(pvalue <= 0.05))

    expect_identical(res_freq@matrix, exp_freq@matrix)
    expect_identical(methods::slotNames(res_freq),
                     methods::slotNames(exp_freq))
})


# Test chord
test_that("Test liana heatmaps", {
    res_chord <- chord_freq(liana_cpdb %>% filter(pvalue <= 0.05))
    expect_identical(class(res_chord), "recordedplot")
})
