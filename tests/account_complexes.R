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

)


# score + columns that it takes into account
# seperate tibble for each score (produced from lr_res)
#


recomplexify <- function(lr_cmplx,
                         entity,
                         columns,
                         funx){

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
                # group_by(source, target, ligand, receptor_complex) %>%
                mutate( {{ entity.col2 }} := min(!!entity.col))

        })

    return(lr_cmplx)
}


recomplexify_check <- function(lr_cmplx,
                               entity,
                               columns,
                               funx){

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
            # group_by(source, target, ligand, receptor_complex) %>%
            mutate( {{ entity.col2 }} := min(!!entity.col))

        })

    return(lr_cmplx)
}





