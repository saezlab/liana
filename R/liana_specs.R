#' S4 Class used to generate aggregate/consesus scores for the methods.
#'
#' @name ScoreSpecifics-class
#'
#' @field method_name name of the method (e.g. cellchat)
#' @field method_score The interaction score provided by the method (typically
#' the score that reflects the specificity of interaction)
#' @field descending_order whether the score should be interpreted in
#'  descending order (i.e. highest score for an interaction is most likely)
#'
#' @field score_fun name of the function to call to generate the results
#'
#' @field columns columns required to generate the score
#'
#' @exportClass ScoreSpecifics
setClass("ScoreSpecifics",
         slots=list(method_name = "character",
                    method_score = "character",
                    descending_order= "logical",
                    score_fun = "function",
                    columns = "character")
)



#' Score Specs Holder
#'
#' @return list of RankSpecifics objects for each method
#'
#' @noRd
#'
#' @details to be explained better and to replace .rank_specs in liana_aggregate
.score_specs <- function(){
    list(
        "connectome" =
            methods::new(
                "ScoreSpecifics",
                method_name = "connectome",
                method_score = "weight_sc",
                descending_order = TRUE,
                score_fun = connectome_score,
                columns = c("ligand.scaled", "receptor.scaled")
            ),
        "logfc" =
            methods::new(
                "ScoreSpecifics",
                method_name = "logfc", # ~italk
                method_score = "logfc_comb",
                descending_order = TRUE,
                score_fun = logfc_score,
                columns = c("ligand.log2FC", "receptor.log2FC")
            ),
        "natmi" =
            methods::new(
                "ScoreSpecifics",
                method_name = "natmi",
                method_score = "edge_specificity",
                descending_order = TRUE,
                score_fun = natmi_score,
                columns = c("ligand.expr", "receptor.expr",
                            "ligand.sum", "receptor.sum")
            )
    )
}


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

#' Rank Specs Holder
#'
#' @return list of RankSpecifics objects for each method
#'
#' @noRd
.rank_specs <- function(){
    list(
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
}