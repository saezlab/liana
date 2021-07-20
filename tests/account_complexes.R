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
#     get_correlation(sce) %>%
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
lr_cdbd %>% distinct_at(vars(ligand, receptor, source,target))



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
           ligand.expr, receptor.expr, ligand.scaled, receptor.scaled) %>%
    # receptor
    group_by(ligand, receptor_complex) %>%
    mutate(receptor.expr.cmplx = min(receptor.expr))# %>%
    top_n(1, receptor.expr.cmplx)
    # ungroup() %>%
    # distinct_at(.vars=c("ligand_complex", "receptor_complex",
                #         "source", "target"),
                # .keep_all=TRUE)


xx2 <- recomplexify_check(lr_cmplx,
                          entity = "receptor",
                          columns = c("receptor.expr",
                                      "receptor.scaled",
                                      "receptor.log2FC",
                                      "ligand.expr",
                                      "ligand.scaled",
                                      "ligand.log2FC"))


# recomplexify should be applied on a per-score basis
# i.e. for NATMI min/mean expr
# min/mean z-score for Conn
# min/mean logFC

score_list <- list(
    edge_specificity = c("receptor.expr", "ligand.expr"),
    weight_sc = c("receptor.scaled", "ligand.scaled"),
    logfc_comb = c("ligand.log2FC", "receptor.log2FC")
)

setClass("ScoreSpecifics",
         slots=list(method_name="character",
                    method_score="character",
                    descending_order="logical",
                    default_fun="function",
                    columns = "character"))



list(
    "connectome" =
        methods::new(
            "ScoreSpecifics",
            method_name = "connectome",
            method_score = "weight_sc",
            descending_order = TRUE,
            default_fun = "",
            columns = c("receptor.scaled", "ligand.scaled"),
        ),
    "logfc_comb" =
        methods::new(
            "ScoreSpecifics",
            method_name = "italk",
            method_score = "logfc_comb",
            descending_order = TRUE,
            default_fun = "",
            columns = c("ligand.log2FC", "receptor.log2FC")
        ),
    "natmi" =
        methods::new(
            "ScoreSpecifics",
            method_name = "natmi",
            method_score = "edge_specificity",
            descending_order = TRUE,
            default_fun = "",
            columns = c("receptor.expr", "ligand.expr"),
        ),
)







connectome_score <- function(lr_res){
    lr_res %>%
        rowwise() %>%
        mutate(weight_sc = mean(c(ligand.scaled, receptor.scaled))) %>%
        select(source, starts_with("ligand"),
               target, starts_with("receptor"),
               everything())

}


#' score + columns that it takes into account
#' seperate tibble for each score (produced from lr_res)



#' Helper function used to account for complexes
#'
#' @param lr_cmplx decomplexified lr_res
#' @param entity ligand or receptor
#' @param columns columns to account for complexes for
#'
#' @returns complex-accounted lr_res
#'
#' @details to be passed before the relevant score_calc function
recomplexify <- function(lr_cmplx,
                         entity,
                         columns){

    # fun = abs+min or mean
    # return the subunit /w min z (by default)

    columns %>%
        map(function(col){
            entity.col = sym(str_glue("{col}"))
            entity.col2 = sym(str_glue("{col}.cmplx")) # to delete

            entity.complex = sym(str_glue("{entity}_complex"))

            alt_entity <- sym(if_else(str_detect(entity, "ligand"),
                                      "receptor", "ligand"))

            lr_cmplx <<- lr_cmplx %>%
                group_by(source, target, !!alt_entity, !!entity.complex) %>%
                mutate( {{ entity.col2 }} := min(!!entity.col))
        })

    return(lr_cmplx)
}


