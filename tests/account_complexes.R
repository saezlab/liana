# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
op_resource <- select_resource("CellPhoneDB")[[1]]

require(SingleCellExperiment)


lr_res <- liana_pipe(seurat_object,
                     op_resource)

xd <- lr_res %>% recomplexify(
    columns = .score_specs()[["logfc"]]@columns,
    complex_policy = "min0")

xd2 <- lr_res %>% recomplexify(
    columns = .score_specs()[["logfc"]]@columns,
    complex_policy = "min")

all_equal(xd, xd2)

xd %>%
    anti_join(xd2, by = c("source", "target",
                          "ligand.complex", "ligand",
                          "receptor.complex", "receptor"))


# if any subunit is missing then assign 0

# take complexes only
# check if a subunit of a complex is missing (rowwise)
# add and assign 0s or pval of 1 to @columns
# then recomplexify


lr_cmplx <- lr_res %>%
    select(ends_with("complex")) %>%
    filter(if_any(.cols = ends_with("complex"),
                  .fns = ~ str_detect(.x, "_")))

complex <- "CNTFR_IL6ST_LIFR"
example <- lr_cmplx %>%
    filter(receptor.complex==complex)

complex_split <- str_split(complex, "_") %>% pluck(1)

absent_subunits <- setdiff(complex_split,
                           lr_res$receptor %>% unique())

# if there are absent subunits assign 0s and pvalues of 1
if(length(absent_subunits)>0){
    lr_res %>%
        # any numeric value to 0
        mutate(across(where(is.numeric) & starts_with("receptor"),
                      ~ ifelse(receptor.complex==complex, 0, .))) %>%
        # FDR and pval to 1
        mutate(across(starts_with("receptor") &
                          (ends_with("FDR") | ends_with("pval")),
                      ~ ifelse(receptor.complex==complex, 1, .)))
}



# Loop over each complex and in lr_res and do the above
lr_cmplx <- lr_res %>%
    select(ends_with("complex")) %>%
    filter(if_any(.cols = ends_with("complex"),
                  .fns = ~ str_detect(.x, "_")))


entity <- "receptor"
entity.complex <- str_glue("{entity}.complex")

complex_entities <- lr_cmplx %>%
    filter(str_detect(.data[[entity.complex]], "_")) %>%
    pluck(entity.complex) %>%
    unique()


# Check if function works
complex <- "CNTFR_IL6ST_LIFR"
example <- lr_res %>%
    filter(receptor.complex==complex)

example_corrected <- example %>%
    missing_subunits_to0(complex = complex,
                         entity = "receptor")
all_equal(example, example_corrected)


# check if it would correct if all units are present
complex <- "ACVR1_TGFBR2"
example <- lr_res %>%
    filter(receptor.complex==complex)

example_corrected <- example %>%
    missing_subunits_to0(complex = complex,
                         entity = "receptor")

all_equal(example, example_corrected)

#
lr_res <- lr_res %>%
    filter(ligand!="ITGAX")





#
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
op_resource <- select_resource("CellPhoneDB")[[1]]

lr_res <- liana_pipe(seurat_object,
                     op_resource)
lr_res2 <- lr_res



recomplexify_env = new.env()
lr_res %<>% account_missing(recomplexify_env)

all_equal(lr_res, lr_res2)




recomplexify()


# Turn into function
account_missing <- function(lr_res, env){

    env$lr_res <- lr_res

    ligand_complexes <- lr_res %>%
        filter(str_detect(.data[["ligand.complex"]], "_")) %>%
        pluck("ligand.complex") %>%
        unique()

    receptor_complexes <- lr_res %>%
        filter(str_detect(.data$receptor.complex, "_")) %>%
        pluck("receptor.complex") %>%
        unique()

    map(receptor_complexes,
        ~missing_subunits_to0(lr_res = lr_res,
                              complex = .x,
                              entity = "receptor",
                              env = env))

    map(receptor_complexes,
        ~missing_subunits_to0(lr_res = lr_res,
                              complex = .x,
                              entity = "ligand",
                              env = env))

    return(env$lr_res)
}






#' Helper Function that assigns 0s to any complexes with missing subunits
#'
#' @param lr_cmplx LR results as obtained by `lr_pipe`
#' @param complex complex of interest
#' @param entity is the complex a 'ligand' or 'receptor'
#'
#' @return A `lr_res` tibble with
missing_subunits_to0 <- function(lr_res, complex, entity, env){

    entity.complex <- str_glue("{entity}.complex")

    complex_split <- # split complex into subunits
        str_split(complex, "_") %>%
        pluck(1)

    # check if subunits are present
    absent_subunits <- setdiff(complex_split,
                               lr_res[[entity]] %>% unique())

    # if there are absent subunits assign 0s and pvalues of 1
    if(length(absent_subunits)>0){
        env$lr_res <- env$lr_res %>%
            # any numeric value to 0
            mutate(across(where(is.numeric) & starts_with(!!entity),
                          ~ ifelse(.data[[entity.complex]]==complex, 0, .))) %>%
            # FDR and pval to 1
            mutate(across(starts_with(!!entity) &
                              (ends_with("FDR") | ends_with("pval")),
                          ~ ifelse(.data[[entity.complex]]==complex, 1, .)))
    }

    return()
}








#
#
#
#
#
#
#
#











# Run /w CellPhoneDB
lr_cdbd <- liana_pipe(seurat_object,
                      select_resource("CellPhoneDB")[[1]])

lr_cdbd



# create a list element with each score
# + custom select for each method
lr_cmplx




scores <- .score_specs() %>%
    map(function(score_object){

        lr_cmplx %<>%  # decomplexify or not
            select(ligand, receptor,
                   ends_with("complex"),
                   source, target,
                   !!score_object@columns) %>%
            recomplexify(columns = score_object@columns) # all or min, ...

        args <- list(
            lr_res = lr_cmplx,
            score_col = score_object@method_score
        )

        exec(score_object@score_fun, !!!args)
})



lr_cmplx %>%
    select(-ends_with("pval")) %>%
    select(-ends_with("FDR")) %>%
    select(-ends_with("stat")) %>%
    select(-ends_with("complex")) %>%
    filter((receptor == "ACVR1" | receptor == "TGFBR2") &
               ligand == "TGFB1") %>%
    filter(source == "B" & target == "B") %>%
    distinct()


# natmi
scores$natmi %>%
    filter(receptor.complex == "ACVR1_TGFBR2" &&
               ligand.complex == "TGFB1") %>%
    filter(source == "B" && target == "B")

# Conn
scores$connectome %>%
    filter(receptor.complex == "ACVR1_TGFBR2" &&
               ligand.complex == "TGFB1") %>%
    filter(source == "B" & target == "B")

# logfc
scores$logfc_comb %>%
    filter(receptor.complex == "ACVR1_TGFBR2" &
               ligand.complex == "TGFB1") %>%
    filter(source == "B" && target == "B")





x2

# scale fun - instead of Seurat?
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)





xd <- recomplexify(lr_cmplx,
                   columns = c("receptor.expr",
                               "receptor.scaled",
                               "receptor.log2FC",
                               "ligand.expr",
                               "ligand.scaled",
                               "ligand.log2FC"))



res1 <- liana_wrap(seurat_object,
                   method = c('connectome',"call_connectome"),
                   resource = c('ICELLNET'),
                   liana_call.params=list("protein"="complex"),
                   decomplexify=TRUE)

res1 <- get_connectome(seurat_object,
                       op_resource = select_resource("ICELLNET")[[1]],
                       decomplexify=TRUE,
                       protein="complex",
                       complex_policy = "min0")

res <- get_connectome(seurat_object,
                      op_resource = select_resource("ICELLNET")[[1]],
                      decomplexify=TRUE,
                      protein="complex",
                      complex_policy = "min0")

res2 <- get_connectome(seurat_object,
                       op_resource = select_resource("ICELLNET")[[1]],
                       decomplexify=TRUE,
                       protein="subunit",
                       complex_policy = "min0")


res1 %>%
    liana_aggregate()

columns


lr_xx <- lr_res %>%
    group_by(across(all_of(c("source", "target",
                             "ligand.complex", "receptor.complex")))) %>%
    # filter(source =="B") %>%
    # filter(target == "B") %>%
    # filter(ligand == "TGFB1") %>%
    select(ligand,receptor,ligand.complex,
           receptor.complex,ligand.scaled,
           receptor.scaled)

xx <- lr_xx %>%
    mutate(receptor.min = min0(receptor.scaled)) %>%
    mutate(receptor.flag = ifelse(receptor.min==receptor.scaled,
                                  TRUE,
                                  FALSE)) %>%
    mutate(ligand.min = min0(receptor.scaled)) %>%
    mutate(ligand.flag = ifelse(receptor.min==receptor.scaled,
                                TRUE,
                                FALSE)) #%>%
# filter(receptor.flag==1 & ligand.flag==1) %>%
# distinct() %>%
# filter(across(ends_with("flag")))


columns %>%
    map(function(col){

        col.min <- sym(str_glue("{col}.min"))
        col.flag <- sym(str_glue("{col}.flag"))

        lr_xx <<- lr_xx %>%
            filter(source =="B") %>%
            filter(target == "B") %>%
            filter(ligand == "TGFB1") %>%
            group_by(across(all_of(c("source", "target",
                                     "ligand.complex",
                                     "receptor.complex")))) %>%
            mutate( {{ col.min }} := min0(.data[[col]])) %>%
            mutate( {{ col.flag }} := ifelse(.data[[col]]==.data[[col.min]],
                                             1,
                                             0))

    })



# returns ties
xx1 <- res %>%
    unite(ligand, receptor,
          source, target, col="xx") %>%
    pluck("xx")
xx2 <- res2 %>% unite(ligand,
                      receptor,
                      source,
                      target,
                      col="xx") %>%
    pluck("xx")
xx3 <- res2 %>% unite(ligand.complex,
                      receptor.complex,
                      source,
                      target,
                      col="xx") %>%
    pluck("xx")


xx3[!(xx3 %in% xx1)]


length(xx2[duplicated(xx2)])



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
lr_cdbd <- liana_pipe(seurat_object,
                      select_resource("CellPhoneDB")[[1]])

lr_cdbd


subunit <- get_connectome(lr_cdbd,
                          protein='subunit')

complex <- get_connectome(lr_cdbd,
                          protein='complex')

complex2 <- get_connectome(lr_cdbd,
                           protein='complex')






# Account for missing subunits
cpdb <- select_resource("CellPhoneDB")[[1]]

cpdb_decomplex <- decomplexify(cpdb)


test <- liana_pipe(seurat_object = testdata,
                   op_resource = cpdb_decomplex)


# xx
cmplx <- op_resource %>%
    select(
        ligand = source_genesymbol,
        ligand.complex = source_genesymbol_complex,
        receptor = target_genesymbol,
        receptor.complex = target_genesymbol_complex
    ) %>%
    filter(if_any(.cols = ends_with("complex"),
                  .fns = ~ str_detect(.x, "_")))

cmplx_og <- op_resource %>%
    select(
        ligand = source_genesymbol,
        ligand.complex = source_genesymbol_complex,
        receptor = target_genesymbol,
        receptor.complex = target_genesymbol_complex
    )

xx <- lr_res %>%
    left_join(., cmplx,
              by=c("ligand", "receptor")) %>%
    distinct()


xx_og <- lr_res %>%
    left_join(., cmplx_og,
              by=c("ligand", "receptor")) %>%
    distinct()
