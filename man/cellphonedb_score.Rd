% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm_funs.R
\name{cellphonedb_score}
\alias{cellphonedb_score}
\title{Function to calculate p-values as in CellPhoneDB}
\usage{
cellphonedb_score(
  lr_res,
  perm_means,
  parallelize,
  workers,
  score_col = "pvalue",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{lr_res}{liana pipe results}

\item{perm_means}{permutations obtained via `get_permutations`}

\item{parallelize}{logical whether to parallelize}

\item{workers}{Number of workers to be used in parallelization}

\item{score_col}{name of the score column}

\item{verbose}{logical for verbosity}

\item{...}{placeholder}
}
\value{
lr_res + pvalue and lr.mean
}
\description{
Function to calculate p-values as in CellPhoneDB
}
\keyword{internal}
