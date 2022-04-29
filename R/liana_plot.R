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
