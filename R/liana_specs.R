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
#' @return list of ScoreSpecifics objects for each method
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
            ),
        "sca" =
            methods::new(
                "ScoreSpecifics",
                method_name = "sca",
                method_score = "LRscore",
                descending_order = TRUE,
                score_fun = sca_score,
                columns = c("ligand.expr", "receptor.expr", "global_mean")
            ),
        "cpdb" =
            methods::new(
                "ScoreSpecifics",
                method_name = "cellphonedb",
                method_score = "pvalue",
                descending_order = FALSE,
                score_fun = cpdb_score,
                columns = c("ligand.trunc", "receptor.trunc")
            ),
        "squidpy" =
            methods::new(
                "ScoreSpecifics",
                method_name = "Squidpy",
                method_score = "pvalue",
                descending_order = FALSE,
                score_fun = function(){},
                columns = ""
            ),
        "cellchat" =
            methods::new(
                "ScoreSpecifics",
                method_name = "cellchat",
                method_score = "pval",
                descending_order = FALSE,
                score_fun = function(){},
                columns = ""
            ),

        # deprecated
        "call_connectome" =
            methods::new(
                "ScoreSpecifics",
                method_name = "connectome",
                method_score = "weight_sc",
                descending_order = TRUE,
                score_fun = function(){},
                columns = ""
            ),
        "call_sca" = methods::new(
            "ScoreSpecifics",
            method_name = "sca",
            method_score = "LRscore",
            descending_order = TRUE,
            score_fun = function(){},
            columns = ""
        ),
        "call_italk" =
            methods::new(
                "ScoreSpecifics",
                method_name = "italk",
                method_score = "logfc_comb",
                descending_order = TRUE,
                score_fun = function(){},
                columns = ""
            ),
        "call_natmi" =
            methods::new(
                "ScoreSpecifics",
                method_name = "natmi",
                method_score = "edge_specificity",
                descending_order = TRUE,
                score_fun = function(){},
                columns = ""
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
                    descending_order="logical")
)
