ligrec <- compile_ligrec_descr()
omni_resources <- compile_ligrec()

CellChatDB <- ligrec[["CellChatDB"]]
CellPhoneDB <- ligrec[["CellPhoneDB"]]

ligrec_decomplexify <- function(ligrec){
    complex_resources <- c("CellChatDB",
                           "CellPhoneDB",
                           "Baccin2019", #*
                           "ICELLNET" #* Complexes present? Why?
                           # CellTalkDB - no complexes in our version
                           )

    cats <- list(transmitters = "genesymbol",
                 receivers = "genesymbol",
                 interactions = c("source_genesymbol",
                                  "target_genesymbol"))

    ligrec <-
        ligrec %>% map2(names(.),
                        function(res, resname) map2(cats, names(cats),
                                                    function(col, cat)
                                                        if(resname %in% complex_resources){
                                                            decomplexify(res[[cat]], col)
                                                            } else{ res[[cat]] }
                                                    )
                        )
    return(ligrec)
}

decomplexify <- function(resource, column){
    column %>%
        map(function(col){
            sep_cols <- c(str_glue("col{rep(1:5)}"))
            resource <<- resource %>%
                separate(col, into = sep_cols, sep = "_") %>%
                pivot_longer(cols = sep_cols,
                             values_to = col,
                             names_to = NULL) %>%
                tidyr::drop_na(col) %>%
                distinct()
        })
    return(resource)
}




