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
                                     "ligand.complex", "receptor.complex")))) %>%
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

