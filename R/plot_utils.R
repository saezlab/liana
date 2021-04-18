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
  }) %>% reduce(., full_join, by = "interaction") %>%
    mutate_at(vars(1:ncol(.)), ~ replace(., is.na(.), 0)) %>%
    mutate_at(vars(2:ncol(.)), ~ replace(., . != 0, 1)) %>%
    as.data.frame()
}


#' Helper function to save Upset plot
#' @param upset_df prepForUpset output
#' @param dir path to save upset plot
#' @return Null
#' @import UpSetR
plotSaveUset <- function(upset_df, dir){
  up <- upset(upset_df, nsets = ncol(upset_df), order.by = "freq",
        point.size = 7, line.size = 2, text.scale	= 2,
        mainbar.y.label = "Significant Interactions Intersect",
        sets.x.label = "Interactions per tool")
  file_name = paste(str_glue(dir))
  png(file_name, width = 1400, height = 900)
  print(up)
  dev.off()
}


#' ??? heatmap
#'
#' @param sig_list named list of significant hits. Named list of methods with
#'     each element being a named list of resources
#' @inheritDotParams pheatmap::pheatmap
#'
#' @return A pheatmap showing binary overlap between methods and resources
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom stringr str_glue
get_BigHeat <- function(sig_list,
                        ...){

  # remove OmniPath and Random resources, as they are much larger than the
  # rest of the resources and result in too much sparsity to get meaningful
  # clusters not completely align to them.
  heatmap_sig_list <- sig_list # %>%
    # map(function(x) x %>%
    #       purrr::list_modify("OmniPath" = NULL) %>%
    #       purrr::list_modify("Random" = NULL)
    # )

  # get method and resource names combined
  lnames <- map(names(heatmap_sig_list), function(m_name){
    map(names(heatmap_sig_list[[m_name]]), function(r_name){
      str_glue("{m_name}_{r_name}")
    })
  }) %>%
    unlist()


  # get binarized significant hits list (1 for sig per method, 0 if absent)
  heatmap_binary_df <- heatmap_sig_list %>%
    purrr::flatten() %>%
    setNames(lnames) %>%
    prepForUpset() %>%
    as_tibble() %>%
    column_to_rownames("interaction")


  # annotation groups (sequential vectors as in heatmap_binary_df)
  method_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "_") %>%
    pull(method)
  resource_groups <- colnames(heatmap_binary_df) %>%
    enframe() %>%
    separate(value, into = c("method", "resource"), sep = "_") %>%
    pull(resource)

  # data frame with column annotations.
  # with a column for resources and a column for methods
  annotations_df <- data.frame(Resource = resource_groups,
                               Method = method_groups)  %>%
    mutate(rn = colnames(heatmap_binary_df)) %>%
    column_to_rownames("rn")

  # List with colors for each annotation.
  mycolors <- list(Method = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(method_groups))),
                   Resource = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(resource_groups))))
  names(mycolors$Resource) <- unique(resource_groups)
  names(mycolors$Method) <- unique(method_groups)

  binary_heatmap <- pheatmap(heatmap_binary_df,
                             annotation_col = annotations_df,
                             annotation_colors = mycolors,
                             ...
                             )

  return(binary_heatmap)
}


#' Helper function to Swap Nested Lists
#' @param sig_list named list of significant hits. Named list of methods with
#' each element being a named list of resources
#' @return A list of resources with each element being a named list of methods
#'
#' @details Swap from nested resource lists to nested method lists,
#'  previously the resources were nested by method, this returns the opposite
get_swapped_list <- function(sig_list){

  sig_df <- sig_list %>%
    enframe() %>%
    unnest(value) %>%
    mutate(name = map(names(sig_list), # get combined method and resource names
                      function(m_name){
                        map(names(sig_list[[m_name]]),
                            function(r_name){
                              str_glue("{m_name}_{r_name}")
                            })
                      }) %>% unlist()) %>%
    separate(name, into = c("method", "resource"), sep = "_") %>%
    mutate(value = value %>% setNames(method)) %>%
    group_by(resource)

  # Keep names of resources
  sig_resource_names <- group_keys(sig_df) %>%
    pull(resource)

  sig_list_resource <- sig_df %>%
    group_split() %>%
    map(function(r_list) # get only list values from resource lists
      r_list %>%
        pull(value)) %>%
    setNames(sig_resource_names)
}




#' PCA plot for Cell type/cluster pair frequency for 'significant/top' hits
#' @param freq_df named list of significant hits. Named list of methods with
#' each element being a named list of resources
#' @return A ggplot2 object
#' @import ggfortify ggplot2 RColorBrewer
plot_freq_pca <- function(freq_df){
  # format to df with frequencies and
  # Resource and Method columns as factors
  cell_pair_frequency <- freq_df %>%
    pivot_wider(names_from = clust_pair,
                values_from = freq,
                id_cols = name,
                values_fill = 0) %>%
    as.data.frame() %>%
    separate(name, into = c("Method", "Resource"), remove = FALSE, sep="_") %>%
    mutate(Method = factor(Method, # prevent ggplot2 from rearranging
                           )) %>%
    mutate(Resource = factor(Resource)) %>%
    column_to_rownames("name")


  # get PCs
  pca_res <- prcomp(cell_pair_frequency[3:ncol(cell_pair_frequency)])

  # frequency plot
  pca_freq <- autoplot(pca_res, data = cell_pair_frequency,
                       colour = "Method", shape = "Resource",
                       size = 6, position = "jitter") +
    scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nlevels(cell_pair_frequency$Method))) +
    theme_bw(base_size = 26) +
    scale_shape_manual(values=1:nlevels(cell_pair_frequency$Resource))

  return(pca_freq)
}

