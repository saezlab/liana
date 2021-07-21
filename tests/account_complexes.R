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


# calculate scores
# lr_res %>%
#     mutate(logfc_comb = product(ligand.log2FC, receptor.log2FC)) %>%
#     get_correlation(sce) %>% # need to change it to be cell-specific
#     # natmi scores
#     rowwise() %>%
#     mutate(edge_specificity = ((ligand.expr*(ligand.sum^-1))) *
#                ((receptor.expr*(receptor.sum^-1)))) %>%
#     mutate(edge_specificity = tidyr::replace_na(edge_specificity, 0)) %>%
#     rowwise() %>%
#     mutate(weight_sc = mean(c(ligand.scaled, receptor.scaled))) %>%
#     select(source, starts_with("ligand"),
#            target, starts_with("receptor"),
#            everything())




# Run /w CellPhoneDB
cpdb <- select_resource("CellPhoneDB")[[1]]
cpdb_decomplex <- cpdb %>%
    decomplexify(column = c("source_genesymbol",
                            "target_genesymbol"))

# Run /w decomplexified cpdb
lr_cdbd <- liana_pipe(seurat_object, cpdb_decomplex)

# Join complexes (recomplexify) to lr_res
cmplx <- cpdb_decomplex %>%
    select(
        ligand = source_genesymbol,
        ligand.complex = source_genesymbol_complex,
        receptor = target_genesymbol,
        receptor.complex = target_genesymbol_complex
    )

lr_cmplx <- lr_cdbd %>%
    left_join(., cmplx,
              by=c("ligand", "receptor")) %>%
    distinct()


# create a list element with each score
# + custom select for each method
lr_cmplx


scores <- .score_specs() %>%
    map(function(score_object){

        args <- list(lr_res = lr_cmplx,
                     score_col = score_object@method_score
                     )

        exec(score_object@score_fun, !!!args) %>%
            select(ligand, receptor,
                   ends_with("complex"),
                   source, target,
                   !!score_object@columns) %>%
            recomplexify(columns = score_object@columns)

    })



scores$connectome



xd <- recomplexify(lr_cmplx,
                   columns = c("receptor.expr",
                               "receptor.scaled",
                               "receptor.log2FC",
                               "ligand.expr",
                               "ligand.scaled",
                               "ligand.log2FC"))



#' Helper function to account for complexes in the resources
#'
#' @param lr_cmplx decomplexified lr_res
#' @param columns columns to account for complexes for
#' @param complex_policy policy how to account for the presence of complexes.
#'   Following the example of \url{https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.ligrec.html}{Squidpy}, valid options are:
#'   'min' select the subunit with the change/expression closest to 0 (as in CellPhoneDB)
#'   'all' returns all subunits seperately
#'
#' @returns complex-accounted lr_res
#'
#' @details to be passed before the relevant score_calc function
#' @importFrom stringr str_split
recomplexify <- function(lr_cmplx,
                         columns){

    # fun = abs+min or mean
    # return the subunit /w min z (by default)

    columns %>%
        map(function(col){

            entity <- as.character(str_split({col}, "\\.")[[1]][[1]])
            entity.col = sym(str_glue("{col}"))
            entity.col2 = sym(str_glue("{col}.cmplx")) # to delete

            entity.complex = sym(str_glue("{entity}.complex"))

            alt_entity <- sym(if_else(str_detect(entity, "ligand"),
                                      "receptor", "ligand"))

            lr_cmplx <<- lr_cmplx %>%
                group_by(source, target, !!alt_entity, !!entity.complex) %>%
                mutate( {{ entity.col2 }} := min(!!entity.col))
        })

    return(lr_cmplx)
}


# lr_acc_cmpx <- lr_cmplx %>%
#     # filter(str_detect(receptor_complex, "_")) %>%
#     select(source, target, ligand, receptor, ligand_complex, receptor_complex,
#            ligand.expr, receptor.expr, ligand.scaled, receptor.scaled) %>%
#     # receptor
#     group_by(ligand, receptor_complex) %>%
#     mutate(receptor.expr.cmplx = min(receptor.expr))# %>%
#     top_n(1, receptor.expr.cmplx)
#     # ungroup() %>%
#     # distinct_at(.vars=c("ligand_complex", "receptor_complex",
#                 #         "source", "target"),
#                 # .keep_all=TRUE)

