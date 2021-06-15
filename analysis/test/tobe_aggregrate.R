#
exp1 <- liana_wrap(seurat_object,
                   method = c('italk', 'sca', 'cellchat',
                              'squidpy', 'natmi', 'connectome'),
                   resource = c('OmniPath', "CellChatDB"),
                   cellchat.params = list(
                       nboot = 100,
                       exclude_anns = NULL,
                       thresh = 1,
                       assay = "RNA",
                       .normalize = TRUE,
                       .do_parallel = FALSE,
                       .raw_use = TRUE
                   ))

temp <- exp1 %>%
    liana_aggregate("OmniPath")


rank_enh(temp$italk.weight_comb, TRUE)





# Separate Function
# Add Robust Ranking Score # in suggests with Conditional

library(RobustRankAggreg)

temp2 <- exp1 %>%
    map(function(res) res %>% pluck("OmniPath")) %>%
    map2(names(.), function(res, method_name){
        method_score <- .rank_specs()[[method_name]]@method_score
        desc_order <- .rank_specs()[[method_name]]@descending_order

        .method = sym(as.character(str_glue("{method_name}.{method_score}")))
        .rank_col = sym(as.character(str_glue("{method_name}.rank")))

        res %>%
            top_n(n=if_else(desc_order,
                            cap,
                            -cap),
                  wt=!!sym(method_score)) %>%
            mutate( {{ .rank_col }} := .rank_enh(.data[[method_score]], desc_order)) %>%
            arrange(!!.rank_col) %>%
            rename( {{ .method }} := method_score) %>%
            select(source, ligand, target, receptor, !!.method, !!.rank_col) %>%
            distinct() %>%
            as_tibble() %>%
            arrange("weight_sc") %>%
            unite(c("source", "ligand",
                    "target", "receptor"), col = "interaction") %>%
            pull("interaction")
    })  %>%
    RobustRankAggreg::aggregateRanks(rmat = rankMatrix(glist),
                                     method = 'stuart')

temp2 <- temp2 %>%
    as_tibble() %>%
    rename(aggregate_rank = Score,
           interaction = Name) %>%
    separate(col = "interaction", sep = "_",
             into = c("source", "ligand", "target", "receptor"))

temp %>% left_join(temp2, by = c("source", "ligand", "target", "receptor"))


temp %>% select(source, ligand, target, receptor, ends_with(".rank"))

glist <- list(italk = exp1$italk %>%
                  arrange(weight_comb) %>%
                  unite(c("source", "ligand",
                          "target", "receptor"), col = "interaction") %>%
                  pull("interaction"),
              sca = exp1$sca %>%
                  arrange("LRscore") %>%
                  unite(c("source", "ligand",
                          "target", "receptor"), col = "interaction") %>%
                  pull("interaction"),
              connectome = exp1$connectome %>%
                  arrange("weight_sc") %>%
                  unite(c("source", "ligand",
                          "target", "receptor"), col = "interaction") %>%
                  pull("interaction"))



xd <- RobustRankAggreg::aggregateRanks(rmat = rankMatrix(glist),
                                       method = 'stuart') %>%
    as_tibble()
xd




# keep only interactions which are in the top 1000 results for each tool
# (to simple function)
temp %>%
    filter_at(vars(ends_with(".rank")), all_vars(. < 200))







glist <- list(sample(letters, 4), sample(letters, 10), sample(letters, 12))
glist %<>% setNames(c("a", "b", "c"))
r = rankMatrix(glist)
r


xd <- RobustRankAggreg::aggregateRanks(rmat = r, method = 'RRA')
xd <- RobustRankAggreg::aggregateRanks(rmat = r, method = 'stuart')
xd <- RobustRankAggreg::aggregateRanks(rmat = r, method = 'mean')
RobustRankAggreg::aggregateRanks(rmat = r, method = 'geom.mean')
RobustRankAggreg::aggregateRanks(rmat = r, method = 'min')


temp %>%
    mutate(stuart_rank =
               select(., ends_with(".rank")))


aggregateRanks(select(temp, ""))

?aggregateRank



# I)
# top_N = lowest method res
# NA = N
# Average(All method ranks)

# II)
#
#
#

temp
# Intersect in top X
temp







# Arrange Helper Function
arrange_enh <- function(.data, wt, descending_order){
    wt <- sym(wt)

    if(descending_order){
        arrange(.data, desc(!!wt))
    } else{
        arrange(.data, !!wt)
    }
}






