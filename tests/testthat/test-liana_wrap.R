# Input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

# Test with OmniPath ----
test_that("Test liana wrapper", {
    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_res.RDS"))
    res1 <- liana_wrap(seurat_object,
                       method = c('logfc','natmi', 'connectome'),
                       resource = c('Consensus'))

    expect_equal(res1, exp1)
})

# Test with Default ----
test_that("Test liana wrapper with default resource", {
    if(!exists("test_external")){
        skip("Not testing externals")
    } else{
        if(!test_external) skip("Not testing externals")
    }
    exp2 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_def_res.RDS"))
    skip_on_ci()
    res2 <- liana_wrap(seurat_object,
                       method = c('sca','call_squidpy', "call_sca"),
                       resource = "Default")

    expect_equal(res2, exp2)
})


test_that("Test liana wrapper return all & supp cols", {
    res3 <- liana_wrap(seurat_object,
                       method = c('logfc','natmi', 'connectome'),
                       resource = c('Consensus'),
                       return_all = TRUE,
                       supp_columns = c("ligand.expr", "receptor.expr",
                                        "ligand.pval", "receptor.FDR"))
    row_number <- res3 %>% map_dbl(function(x) nrow(x)) %>% mean()
    expect_equal(4770, row_number)

    row_number_filt <- res3 %>% map_dbl(function(x)
        x %>% filter(!to_filter) %>%
        nrow()) %>%
        mean()
    expect_equal(735, row_number_filt)
})


# Test Default Params ----
test_that("Test liana wrapper", {
    exp3 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_def_args.RDS"))
    res3 <- liana_defaults(expr_prop=0,
                           squidpy.params=list(threshold = 0.1),
                           cellchat.params=list(nboot=1000))

    expect_equal(res3, exp3)
})

# Test Expression_prop filtering ----
test_that("Test expr_prop filtering", {
    liana_pipe_res <- readRDS(file.path(liana_path, "testdata",
                                        "output", "liana_pipe.RDS"))

    expect_equal(nrow(.filt_liana_pipe(liana_pipe_res, "connectome", expr_prop=0)), 5076)
    expect_equal(nrow(.filt_liana_pipe(liana_pipe_res, "cellphonedb")), 996)

})

# Test /w weights/contraints  ----
test_that("Test cell.adj weights", {
    # Generate putative cell adjacency score
    # Here we generate weights between 0 to 1, i.e. similar to densities
    # We assign 0s to cell types which are not expected to communicate
    set.seed(1)
    cell.adj <- readRDS(file.path(liana_path,
                                  "testdata",
                                  "output",
                                  "liana_pipe.RDS")) %>%
        select(source, target) %>%
        distinct() %>%
        mutate(adjacency = rbinom(n=nrow(.), size=1, prob=0.5)) %>%
        rowwise() %>%
        mutate(adjacency = if_else(adjacency==1, runif(1, min=0, max=1), 0)) %>%
        ungroup()

    exp1 <- readRDS(file.path(liana_path, "testdata",
                              "output", "liana_constrained.RDS"))
    res1 <- liana_wrap(seurat_object,
                       method = c('logfc','sca', 'connectome'),
                       resource = c('Consensus'),
                       cell.adj = cell.adj)
    expect_equal(res1, exp1)
})
