% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liana_plot.R
\name{liana_heatmap}
\alias{liana_heatmap}
\title{Communication by cell type Heatmap}
\usage{
liana_heatmap(
  mat,
  font_size = 12,
  grid_text = FALSE,
  name = "Frequency",
  pallette = c("white", "violetred2"),
  row_title = "Sender (Cell types)",
  column_title = "Receiver (Cell types)",
  ...
)
}
\arguments{
\item{mat}{Diagonal celltype-celltype matrix to be plotted. In theory,
any metric deemed meaningful between cell pairs can be plotted.}

\item{font_size}{base font_size - other fontsizes are relative to this one}

\item{grid_text}{logical whether to display grid text or not}

\item{name}{name of the heatmap.
By default the heatmap name is used as the title of the heatmap legend.}

\item{row_title}{Row tittle}

\item{column_title}{Column tittle}

\item{...}{parameters passed to `ComplexHeatmap::Heatmap`}
}
\description{
Communication by cell type Heatmap
}
\details{
Heatmap function inspired by CellPhoneDBv3 and CellChat's designs
on communication heatmaps.
}
