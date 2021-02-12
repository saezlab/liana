#' Helper function to prepare list for UpsetPlot
#' @param named_list a named list with LR results
prepForUpset <- function(named_list){
  map(names(named_list), function(l_name){
      named_list[[l_name]] %>%
      select(1:4) %>%
      unite("interaction", source, target,
            ligand, receptor, sep=".") %>%
      mutate(!!l_name := 1)
  }) %>% reduce(., full_join) %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . > 0, 1)) %>%
    as.data.frame()
}
