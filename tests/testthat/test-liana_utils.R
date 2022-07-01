# Input
liana_path <- system.file(package = "liana")
op_resource <- select_resource('Consensus')[[1]] %>% head(50)

symbols_dict <- readRDS(file.path(liana_path, "human_mouse_orthologues.RDS"))

# Test Liana Pipe
test_that("Test liana pipe", {

    op_ortho <- generate_homologs(op_resource = op_resource,
                                  symbols_dict = symbols_dict)
    res1 <- readRDS(file.path(liana_path,
                              "testdata",
                              "output",
                              "op_ortho.RDS"))

    expect_equal(res1, op_ortho)
})
