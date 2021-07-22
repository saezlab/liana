# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

require(SingleCellExperiment)


res1 <- liana_wrap(seurat_object,
                   method = c('connectome',"call_connectome"),
                   resource = c('CellPhoneDB'),
                   liana_call.params=list("protein"="complex"),
                   decomplexify=TRUE)

res1 <- get_connectome(seurat_object,
               op_resource = select_resource("CellPhoneDB")[[1]],
               decomplexify=TRUE,
               protein="complex",
               complex_policy = "min0")

res <- get_connectome(seurat_object,
                      op_resource = select_resource("CellPhoneDB")[[1]],
                      decomplexify=TRUE,
                      protein="complex",
                      complex_policy = "min0")

res2 <- get_connectome(seurat_object,
                       op_resource = select_resource("CellPhoneDB")[[1]],
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
