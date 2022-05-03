#' Liana dotplot interactions by source and target cells
#'
#' @param liana_agg aggregated `liana_wrap` results -> preferentially filtered
#' by some condition (e.g. preferential ranking, specific interactions, etc)
#'
#' @param source_groups names of the source(sender) cell types
#' @param target_groups names of the target cell types
#'
#' @param magnitude column to represent interactions expression magnitude
#' (by default `sca.LRscore`)
#'
#' @param specificity column to represent the dot-size of the interaction
#' (by default `natmi.edge_specificity`)
#'
#' @param y.label y label name
#' @param size.label size (~specificty) label name
#' @param colour.label colour (~magnitude) label name
#'
#' @param show_complex logical whether to show complexes (default - TRUE) or
#'  only the subunit with minimum expression.
#'
#' @details Here, we refer to `specificity` as how specific this interaction is
#' to a cell type pair regards to the rest of the cell type pairs (
#' e.g. CellPhoneDB's p-values, NATMI's specificity edges, Connectome's scaled weights, etc)
#'
#' `magnitude` on the other hand is a direct measure of the expression alone,
#' by default we use SingleCellSignalR's dataset indepent LRscore (bound between 0 and 1).
#' Yet, one could also use CellChat's probabilities or CellPhoneDB's means, etc.
#'
#' @import ggplot2 dplyr
#'
#' @return a ggplot2 object
#'
#' @export
liana_dotplot <- function(liana_agg,
                          source_groups,
                          target_groups,
                          specificity = "natmi.edge_specificity",
                          magnitude = "sca.LRscore",
                          y.label = "Interactions (Ligand -> Receptor)",
                          size.label = "Interaction\nSpecificity",
                          colour.label = "Expression\nMagnitude",
                          show_complex = TRUE){

    if(show_complex){
        entities <- c("ligand.complex", "receptor.complex")
    } else{
        entities <- c("ligand", "receptor")
    }

    # Modify for the plot
    liana_mod <- liana_agg %>%
        # Filter to only the cells of interest
        filter(source %in% source_groups) %>%
        filter(target %in% target_groups) %>%
        rename(magnitude = !!magnitude) %>%
        rename(specificity = !!specificity) %>%
        unite(entities, col = "interaction", sep = " -> ") %>%
        unite(c("source", "target"), col = "source_target", remove = FALSE)

    # colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    cbPalette <- c("#E69F00", "#56B4E9",
                   "#009E73", "#F0E442", "#0072B2",
                   "#D55E00", "#CC79A7")

    # plot
    suppressWarnings(
        ggplot(liana_mod,
               aes(x = target,
                   y = interaction,
                   colour = magnitude,
                   size = specificity,
                   group = target
               )) +
            geom_point() +
            scale_color_gradientn(colours = viridis::viridis(20)) +
            scale_size_continuous(range = c(5, 9)) +
            facet_grid(. ~ source,
                       space = "free",
                       scales ="free",
                       switch="y") +
            theme_bw(base_size = 20) +
            theme(
                legend.text = element_text(size = 16),
                axis.text.x = element_text(colour =
                                               cbPalette[1:length(
                                                   unique(liana_mod$source)
                                               )],
                                           face = "bold",
                                           size = 23),
                axis.text.y = element_text(size = 18,
                                           vjust = 0.5),
                legend.title = element_text(size = 18),
                panel.spacing = unit(0.1, "lines"),
                strip.background = element_rect(fill = NA),
                strip.text = element_text(size = 24, colour = "gray6") #,
                # strip.text.y.left = element_text(angle = 0)
            ) +
            # scale_x_discrete(position = "right") +
            labs(y = y.label,
                 colour = colour.label,
                 size = size.label,
                 x = NULL
            )
    )
}




#' Communication Frequency heatmap plot
#'
#' @param liana_res aggregated liana results (preferably truncated
#'  to some threshold)
#' @param pallette colour palette to be used for the cells
#' @inheritDotParams liana_heatmap
#'
#' @export
#'
#' @details This plot was inspired by CellPhoneDB and also CellChat's heatmap design.
#' It makes the assumption that the number of interactions inferred between cell
#' types is informative of the communication events occurring in the system as a whole.
#' This is a rather strong assumption limited by the arbitrarily filtered
#' interactions Thus, I suggest that one limits any conclusions, unless supported
#' by complimentary information, such as biological prior knowledge.
heat_freq <- function(liana_res, ...){
    # Calculate Frequencies
    freqs <- liana_res %>%
        group_by(source, target) %>%
        summarise(freq = n(), .groups = 'keep') %>%
        pivot_wider(id_cols = source,
                    names_from = target,
                    values_from = freq,
                    values_fill = 0) %>%
        arrange(source) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames('source') %>%
        as.matrix()

    liana_heatmap(mat = freqs,
                  ...)
}


#' Communication Specificity heatmap plot
#'
#' @param liana_res aggregated liana results (preferably truncated
#'  to some threshold)
#' @inheritParams liana_heatmap
#' @param pallette colour palette to be used for the cells
#' @inheritDotParams liana_heatmap
#'
#' @export
#'
#' @details This function builds on the frequency heatmap. Here, we plot the
#' average specificity per celltype pair, using NATMI's specificity edges.
#' NATMI's specificity ranges from 1 to 0, with 1 being assigned to interactions
#' uniquely expressed in a pair of cell types. Similarly to the previous plot,
#' this one is also limited by the arbitrarity of filtering and any conclusions
#' should be considered with caution.
heat_spec <- function(liana_res,
                      name = 'Average\nSpecificity',
                      pallette = c("white", "#2cb7b9"),
                      ...){

    # Average Specificity
    avg_specificity <- liana_trunc %>%
        select(source, target, natmi.edge_specificity) %>%
        group_by(source, target) %>%
        summarise(avg_spec = mean(natmi.edge_specificity),
                  .groups = 'keep') %>%
        arrange(source) %>%
        ungroup() %>%
        pivot_wider(id_cols = source,
                    names_from = target,
                    values_from = avg_spec,
                    values_fill = 0) %>%
        arrange(source) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames('source') %>%
        as.matrix()

    liana_heatmap(mat = avg_specificity,
                  name = name,
                  pallette = pallette,
                  ...)
}





#' Communication by cell type Heatmap
#'
#' @param mat Diagonal celltype-celltype matrix to be plotted. In theory,
#' any metric deemed meaningful between cell pairs can be plotted.
#' @param font_size base font_size - other fontsizes are relative to this one
#' @param grid_text logical whether to display grid text or not
#' @param name name of the heatmap.
#' By default the heatmap name is used as the title of the heatmap legend.
#' @param row_title Row tittle
#' @param column_title Column tittle
#'
#' @param ... parameters passed to `ComplexHeatmap::Heatmap`
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom grid gpar unit grid.text
#'
#' @details Heatmap function inspired by CellPhoneDBv3 and CellChat's designs
#' on communication heatmaps.
liana_heatmap <- function(mat,
                          font_size = 12,
                          grid_text = FALSE,
                          name = 'Frequency',
                          pallette = c("white", "violetred2"),
                          row_title = "Sender (Cell types)",
                          column_title = "Receiver (Cell types)",
                          ...){

    if(grid_text){
        grid_text <- function(j, i, x, y, width, height, fill) {
            grid_text <- grid.text(sprintf("%d", mat[i, j]),
                                   x, y, gp = gpar(fontsize = font_size*0.83))
        }
    } else{
        grid_text <- NULL
    }

    # define Annotations and Barplots
    cell_anno <- unique(rownames(mat))
    cell_anno <- grDevices::colorRampPalette(
        (RColorBrewer::brewer.pal(n = 8, name = 'Dark2'))
    )(length(cell_anno)) %>%
        setNames(cell_anno)

    ## Annotations
    ha_opts <- list(show_legend = FALSE,
                    show_annotation_name = FALSE,
                    col = list("anno"=cell_anno),
                    simple_anno_size = grid::unit(0.25, "cm"))
    column_ha <- exec("HeatmapAnnotation", anno = names(cell_anno), !!!ha_opts)
    row_ha <- exec("rowAnnotation", anno = names(cell_anno), !!!ha_opts)

    # Barplots
    column_bar <- ComplexHeatmap::HeatmapAnnotation(
        bar = .anno_barplot(colSums(mat),
                            cell_anno,
                            axis.font.size = font_size*0.5
                            ),
        annotation_name_gp = gpar(fontsize = font_size*0.5),
        show_legend = FALSE,
        show_annotation_name = FALSE)

    row_bar <- ComplexHeatmap::rowAnnotation(
        bar2 = .anno_barplot(rowSums(mat),
                             cell_anno,
                             font_size*0.5
                             ),
        gp = gpar(fill = cell_anno,
                  col = cell_anno),
        show_legend = FALSE,
        show_annotation_name = FALSE)

    # Heatmap
    ComplexHeatmap::Heatmap(mat,
                            col=colorRampPalette(pallette)(10),
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            row_names_side = "left",
                            top_annotation = column_bar,
                            bottom_annotation = column_ha,
                            right_annotation = row_bar,
                            left_annotation = row_ha,
                            row_title = row_title,
                            row_names_gp = gpar(fontsize = font_size),
                            row_title_gp = gpar(fontsize = font_size*1.2),
                            column_names_gp = gpar(fontsize = font_size),
                            column_title = column_title,
                            column_title_gp = gpar(fontsize = font_size*1.2),
                            column_title_side = "bottom",
                            heatmap_legend_param = list(title_gp = gpar(fontsize = font_size*0.9,
                                                                        fontface = 'bold'),
                                                        border = NA,
                                                        labels_gp = gpar(fontsize = font_size*0.9),
                                                        grid_width = unit(2, "mm")),
                            name = name,
                            cell_fun = grid_text,
                            ...
    )

}


#' Helper Function to Generate Annotation Barplots
#'
#' @param x numeric vector
#' @param cell_anno vector of colour codes named by cell type annotations
#' @param axis.font.size fontsize of the barplots axis font size
#'
#' @noRd
.anno_barplot <- function(x,
                          cell_anno,
                          axis.font.size){
    anno_barplot(x,
                 gp = gpar(fill = cell_anno,
                           col = cell_anno,
                           font.size=axis.font.size),
                 axis_param = list(gp=gpar(font.size=axis.font.size)),
                 title="",
                 border = FALSE)
}



