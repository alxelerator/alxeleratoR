% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permanova2.R
\name{permanova2}
\alias{permanova2}
\title{Single (multivariable) PERMANOVA v2 Uses vegan::adonis2. This function is just a convenience wrapper for vegan::adonis2.}
\usage{
permanova2(
  physeq = "",
  distm = "bray",
  formula = "",
  var = "",
  permutations = 999,
  hellinger = TRUE,
  rngseed = 2202,
  ...
)
}
\arguments{
\item{physeq}{phyloseq object. A distance matrix is created using \code{otu_table(physeq)}; permanovas are computed for each element of \code{sample_variables(physeq)}}

\item{distm}{choice of dissimilarity index according to \code{\link{vegdist}}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".}

\item{formula}{Defaults to distnmatrix ~ var}

\item{var}{ONE variable. Required if no formula is specified! (should be factor levels >=2).}

\item{hellinger}{TRUE/FALSE. By default performs a Hellinger transformation of the OTU table. Here you can turn it off.}

\item{rngseed}{set.seed number to use for reproducible permutations. Default rngseed=2202.}

\item{...}{Allows additional arguments to the adonis2 function to be taken along. For instance parallel = parallel::detectCores() - 1}

\item{perms}{the number of permutations to be performed. Default is 999}
}
\value{
vegan2 default output.
}
\description{
Function performs individual, \code{\link{adonis2}} permanova for given variable or (multivariable) formula
}
\note{
if \code{sample_data(physeq)} contains NAs, the function will convert these to an additional level "none"
}
\examples{
  per <- permanova2( physeq, distm="bray", formula="~ Treatment + AnimalNumber", var="", rngseed="2202", parallel=(parallel::detectCores() - 1) )
  per <- permanova2( physeq, distm="bray", formula="", var="Treatment", rngseed="2202")

}
\author{
Convenience wrapper: Alex Bossers \email{a.bossers@uu.nl}
}
