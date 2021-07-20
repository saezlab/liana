#' Function Used to Calculate the Connectome-like `weight_sc` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#'
#' @return lr_res with an added `weight_sc` column
connectome_score <- function(lr_res){
    lr_res %>%
        rowwise() %>%
        mutate(weight_sc = mean(c(ligand.scaled, receptor.scaled))) %>%
        select(source, starts_with("ligand"),
               target, starts_with("receptor"),
               everything())
}


#' Function Used to Calculate the Connectome-like `weight_sc` weights
#'
#' @param lr_res \link(liana::liana_pipe) results
#'
#' @return lr_res with an added `weight_sc` column
natmi_score <- function(lr_res){
    lr_res %>%
        rowwise() %>%
        mutate(weight_sc = mean(c(ligand.scaled, receptor.scaled))) %>%
        select(source, starts_with("ligand"),
               target, starts_with("receptor"),
               everything())
}









