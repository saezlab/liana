
.rank_cap = 500 # custom N
.rank_cap = exp1 %>% map(function(res) nrow(res)) %>% as.numeric %>% min # optional
.rank_cap = exp1 %>% map(function(res) nrow(res)) %>% as.numeric %>% max # default


# aggregate results
exp1 %>%
    map2(names(.), function(res, method_name){

        if(!is_tibble(res[[1]]) && is.null(resource)){
            stop("Please provide provide a name for the resource, ",
                 "otherwise the first resource will be plucked!")
        } else if(!is_tibble(res[[1]])){
            res %<>% map(function(m_results) m_results %>% pluck(resource))
        }

        method_score <- .rank_specs[[method_name]]@method_score
        desc_order <- .rank_specs[[method_name]]@descending_order
        .method = sym(as.character(str_glue("{method_name}.{method_score}")))
        .rank_col = sym(as.character(str_glue("{method_name}.rank")))

        cap %<>% `%||%`(.select_cap(liana_res, set_cap))

        res %>%
            .list2tib %>%
            top_n(n=if_else(desc_order,
                            .rank_cap,
                            -.rank_cap),
                    wt=!!sym(method_score)) %>%
            mutate({{ .rank_col }} := rank_enh(!!sym(method_score), desc_order)) %>%
            arrange(!!.rank_col) %>%
            rename( {{ .method }} := method_score) %>%
            select(source, ligand, target, receptor, !!.method, !!.rank_col) %>%
            distinct() %>%
            as_tibble()
    }) %>%
    purrr::reduce(., full_join, by = c("source", "ligand", # Full join all results
                                       "target", "receptor")) %>%

    # Better to be a separate function
    mutate_at(vars(ends_with(".rank")),
              ~ replace(., is.na(.), .rank_cap)) %>% # assign .rank_cap to NA
    mutate(min_rank = pmap_dbl(select(., ends_with(".rank")),
                               function(...) min(c(...))),
           mean_rank = pmap_dbl(select(., ends_with(".rank")),
                                function(...) mean(c(...))),
           median_rank = pmap_dbl(select(., ends_with(".rank")),
                                  function(...) median(c(...)))) %>%
    arrange(mean_rank) %>%
    select(source, ligand, target, receptor, ends_with("_rank"), everything())
temp



exp1 %>%
    map2(names(.), function(res, method_name){
        res %>% .list2tib
    })

exp1 <- liana_wrap(seurat_object,
                   method = c('italk', 'sca'),
                   resource = c('OmniPath', "CellChatDB"))

temp <- exp1 %>%
    liana_aggregate("OmniPath")


rank_enh(temp$italk.weight_comb, TRUE)





# keep only interactions which are in the top 1000 results for each tool
# (to simple function)
temp %>%
    filter_at(vars(ends_with(".rank")), all_vars(. < 200))



# Separate Function
# Add Robust Ranking Score # in suggests with Conditional
library(RobustRankAggreg)
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






