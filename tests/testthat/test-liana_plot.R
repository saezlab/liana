# Input
liana_path <- system.file(package = "liana")
liana_cpdb <- readRDS(file.path(liana_path, "testdata",
                                "output", "liana_cpdb.RDS"))

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
