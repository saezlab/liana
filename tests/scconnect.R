liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata",
                      "input", "testdata.rds"))

# liana Pipe Output ----
pipe_out <- liana_pipe(seurat_object,
                       op_resource = select_resource("OmniPath")[[1]] %>%
                           liana:::decomplexify())

scconnect_score(pipe_out %>%
                    filter(ligand.prop >= 0.1) %>%
                    filter(receptor.prop >= 0.1))


pipe_out %>%
    filter(ligand.prop >= 0.1) %>%
    filter(receptor.prop >= 0.1) %>%
    select(ligand, receptor,
           contains("complex"),
           ligand.expr, receptor.expr,
           ligand.pval, receptor.pval) %>%
    rowwise() %>%
    mutate(specificity = -log10(mean(c(ligand.pval, receptor.pval)))) %>%
    mutate(score = sqrt(ligand.expr * receptor.expr)) %>%
    mutate(importance = log10(score) * specificity) %>%
    arrange(desc(importance))


# %>%
    # select(ligand, receptor, contains("complex"), importance)



liana_wrap(seurat_object, method="scconnect")
