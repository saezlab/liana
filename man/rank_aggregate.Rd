% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liana_aggregate.R
\name{rank_aggregate}
\alias{rank_aggregate}
\title{Aggregate CCC Method results and by both magnitude and specificity ranks}
\usage{
rank_aggregate(liana_res, ...)
}
\arguments{
\item{liana_res}{LIANA results}

\item{...}{
  Arguments passed on to \code{\link[=liana_aggregate]{liana_aggregate}}
  \describe{
    \item{\code{aggregate_how}}{way to aggregate, by default (NULL) will aggregate
all passed methods with the approach specified in `liana:::.score_specs`.
Alternative options are `magnitude` and `specificity`.}
    \item{\code{set_cap}}{Function used to set ranked cap (i.e. the value that is
assigned to interactions with NA for scores);
By default, this is set to "max", which is the maximum number of interactions
obtained by the methods; Some methods return all possible ligand-receptor
combinations for each possible source and target cell pair - i.e. the
known universe of all possible interactions (based on the CCC resource)}
    \item{\code{resource}}{If methods are ran with multiple resources, the name of the
resource of interest needs to be provided
*Note* if a name is not provided, the first results based on the first
resource in the list will be returned}
    \item{\code{cap}}{A cap can for all methods can also be manually set, then the top X
interactions, based on the `specificity` scores for each method will be
returned and the ranking will be carried out solely on them}
    \item{\code{get_ranks}}{boolean, whether to return consensus ranks for methods}
    \item{\code{get_agrank}}{boolean, whether to return aggregate rank using the
`RobustRankAggreg` package.}
    \item{\code{.score_mode}}{defines the way that the methods would be aggragate.
By default, we use the score of each method which reflects specificity
(if available), if not e.g. the case of SCA we use it's sole scoring function.
This aggregation is by default done on the basis of the list returns by
`.score_mode`. Alternatively, one could pass `.score_housekeep` to obtain an
aggragate of the housekeeping interactions of each method.}
    \item{\code{join_cols}}{columns by which different method results will be joined.
NULL by default, and automatically will handle the columns depending on the
methods used.}
  }}
}
\description{
Aggregate CCC Method results and by both magnitude and specificity ranks
}
\examples{
liana_path <- system.file(package = "liana")
# load testdata
testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
# run liana
liana_res <- liana_wrap(testdata, method=c("sca", "natmi"))
# aggregate results from multiple methods
liana_res <- rank_aggregate(liana_res)
}
