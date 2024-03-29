% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liana_pipe.R
\name{join_means}
\alias{join_means}
\title{Join Expression per Cluster}
\usage{
join_means(lr_res, means, source_target, entity, type, pb = NULL)
}
\arguments{
\item{lr_res}{LR formatted DE results from \link{ligrec_degformat}}

\item{means}{Gene avg expression per cluster}

\item{source_target}{target or source cell}

\item{entity}{ligand or receptor}

\item{type}{type of mean to join (count or scaled)}

\item{pb}{progress bar}
}
\value{
Returns the Average Expression Per Cluster
}
\description{
Join Expression per Cluster
}
\keyword{internal}
