# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
op_resource <- select_resource("OmniPath")[[1]]

require(SingleCellExperiment)

# Run /w OmniPath
lr_res <- liana_pipe(seurat_object,
                     op_resource)
lr_res


# Run /w CellPhoneDB
cpdb <- select_resource("CellPhoneDB")[[1]]
cpdb_decomplex <- cpdb %>%
    decomplexify(column = c("source_genesymbol",
                            "target_genesymbol"))

# Run /w decomplexified cpdb
lr_cdbd <- liana_pipe(seurat_object, cpdb_decomplex)
lr_cdbd %>% distinct_at(vars(ligand,receptor,source,target))



# recomplexify
cmplx <- cpdb_decomplex %>%
    select(
        ligand = source_genesymbol,
        ligand_complex = source_genesymbol_complex,
        receptor = target_genesymbol,
        receptor_complex = target_genesymbol_complex)

lr_cmplx <- lr_cdbd %>%
    left_join(., cmplx,
              by=c("ligand", "receptor")) %>%
    distinct()


xd <- lr_cmplx %>%
    # filter(str_detect(receptor_complex, "_")) %>%
    select(source, target, ligand, receptor, ligand_complex, receptor_complex,
           ligand.expr, receptor.expr) %>%
    # receptor
    group_by(ligand, receptor_complex) %>%
    mutate(receptor.expr.cmplx = min(receptor.expr)) %>%
    top_n(1, receptor.expr.cmplx) %>%
    ungroup() %>%
    # select(-c(ligand, receptor)) %>%
    # rename(ligand = ligand_complex,
    #        receptor = receptor_complex) %>%
    distinct_at(.vars=c("ligand_complex", "receptor_complex", "source", "target"), .keep_all=TRUE)
    # ligand
    # group_by(ligand_complex, receptor) %>%
    # mutate(ligand.expr.cmplx = min(ligand.expr)) %>%
    # mutate(ligand.rank = rank(ligand.expr.cmplx)) %>%
    # ungroup()


# For each subunit in a complex set the abs minimum/mean of scaled, logFC, expr





#' Helper Function to 'decomplexify' ligands and receptors into
#' @param resource a ligrec resource
#' @param column column to separate and pivot long (e.g. genesymbol or uniprot)
#'
#' @return returns a longer tibble with complex subunits on seperate rows
decomplexify <- function(resource, column){
    column %>%
        map(function(col){
            sep_cols <- c(str_glue("col{rep(1:5)}"))
            col.complex <- str_glue("{col}_complex")

            resource <<- resource %>%
                mutate({{ col.complex }} :=
                           resource[[str_glue("{col}")]]) %>%
                separate(col,
                         into = sep_cols,
                         sep = "_",
                         extra = "drop",
                         fill = "right") %>%
                pivot_longer(cols = sep_cols,
                             values_to = col,
                             names_to = NULL) %>%
                tidyr::drop_na(col) %>%
                distinct() %>%
                mutate_at(.vars = c(col),
                          ~str_replace(., "COMPLEX:", ""))
        })
    return(resource)
}
