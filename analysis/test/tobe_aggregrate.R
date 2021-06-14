
.rank_cap = 500 # custom N
.rank_cap = exp1 %>% map(function(res) nrow(res)) %>% as.numeric %>% min # optional
.rank_cap = exp1 %>% map(function(res) nrow(res)) %>% as.numeric %>% max # default


# aggregate results
temp <- exp1 %>%
    map2(names(.), function(res, method_name){
        # pluck specific resource?
        method_score <- .rank_specs[[method_name]]@method_score
        parm_order <- .rank_specs[[method_name]]@descending_order
        .method = sym(as.character(str_glue("{method_name}.{method_score}")))
        .rank_col = sym(as.character(str_glue("{method_name}.rank")))

        res %>%
            .list2tib %>%
            top_n(n=if_else(parm_order,
                            .rank_cap,
                            -.rank_cap),
                    wt=!!sym(method_score)) %>%
            mutate({{ .rank_col }} := dense_rank(desc(!!sym(method_score)))) %>%
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

is_tibble(exp1[[1]])
is_tibble(exp1)

exp1 %>%
    liana_aggregate("OmniPath")

# keep only interactions which are in the top 1000 results for each tool
# (to simple function)
temp %>%
    filter_at(vars(ends_with(".rank")), all_vars(. < 1000))



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





#' S4 Class used to generate aggregate/consesus scores for the methods.
#'
#' @name RankSpecifics-class
#'
#' @field method_name name of the method (e.g. cellchat)
#' @field method_score The interaction score provided by the method (typically
#' the score that reflects the specificity of interaction)
#' @field descending_order whether the score should be interpreted in
#'  descending order (i.e. highest score for an interaction is most likely)
#'
#' @exportClass RankSpecifics
setClass("RankSpecifics",
         slots=list(method_name="character",
                    method_score="character",
                    descending_order="logical"))

# Rank Specifics Holder
.rank_specs <- list(
    "cellchat" =
        methods::new(
            "RankSpecifics",
            method_name = "cellchat",
            method_score = "pval",
            descending_order = FALSE
        ),
    "connectome" =
        methods::new(
            "RankSpecifics",
            method_name = "connectome",
            method_score = "weight_sc",
            descending_order = TRUE
        ),
    "italk" =
        methods::new(
            "RankSpecifics",
            method_name = "italk",
            method_score = "weight_comb",
            descending_order = TRUE
        ),
    "natmi" =
        methods::new(
            "RankSpecifics",
            method_name = "natmi",
            method_score = "edge_specificity",
            descending_order = TRUE
        ),
    "sca" = methods::new(
        "RankSpecifics",
        method_name = "sca",
        method_score = "LRscore",
        descending_order = TRUE
    ),
    "squidpy" =
        methods::new(
            "RankSpecifics",
            method_name = "Squidpy",
            method_score = "pvalue",
            descending_order = FALSE
        )
)


# Arrange Helper Function
arrange_enh <- function(.data, wt, descending_order){
    wt <- sym(wt)

    if(descending_order){
        arrange(.data, desc(!!wt))
    } else{
        arrange(.data, desc(!!wt))
    }
}






