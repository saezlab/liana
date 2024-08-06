# Input
liana_path <- system.file(package = "liana")
op_resource <- select_resource('Consensus')[[1]] %>% head(50)

sce  <- readRDS(file.path(liana_path , "testdata",
                          "input", "testsce.rds"))
# add random samples
set.seed(1)
sce@colData$sample <- sample(LETTERS[4:5], 90, replace=TRUE)

symbols_dict <- readRDS(file.path(liana_path, "human_mouse_orthologues.RDS"))

# Test Liana Pipe
test_that("Test generate_homologs", {
    op_ortho <- generate_homologs(op_resource = op_resource,
                                  symbols_dict = symbols_dict,
                                  target_organism = "mouse")
    res1 <- readRDS(file.path(liana_path,
                              "testdata",
                              "output",
                              "op_ortho.RDS"))
})


test_that("Test test liana by sample and rank_aggregate", {
    sce <- liana_bysample(sce,
                          idents_col="label",
                          sample_col="sample",
                          permutation.params = list(nperms = 2),
                          aggregate_how = "both"
        )
    expect_equal(c("D", "E"), names(sce@metadata$liana_res))

    expect_true(all(c("magnitude_rank", "specificity_rank") %in%
                    colnames(sce@metadata$liana_res$D)))

    best_inverted <- sce@metadata$liana_res %>%
        preprocess_scores(score_col = "magnitude_rank") %>%
        pluck("E") %>%
        pluck("magnitude_rank") %>%
        pluck(1)
    expect_equal(best_inverted, 1)
})
