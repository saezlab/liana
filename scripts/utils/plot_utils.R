#' Helper function to prepare list for UpsetPlot
#' @param named_list a named list with sugnificant LR results
#' @return Matrix/DF of 0 and 1 where 1 means that the interaction
#' is present and 0 means that it is not
prepForUpset <- function(named_list){
  map(names(named_list), function(l_name){
      named_list[[l_name]] %>%
      select(1:4) %>%
      unite("interaction", source, target,
            ligand, receptor, sep="_") %>%
      mutate(!!l_name := 1)
  }) %>% reduce(., full_join) %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
    as.data.frame()
}



plotSaveUset <- function(upset_df, dir){
  up <- upset(upset_df, nsets = ncol(upset_df), order.by = "freq",
        point.size = 5, line.size = 2, text.scale	= 2,
        mainbar.y.label = "Significant Interactions Intersect",
        sets.x.label = "Interactions per tool")
  file_name = paste(str_glue(dir))
  png(file_name, width = 1200, height = 900)
  print(up)
  dev.off()
}
