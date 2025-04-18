% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permanovas.R
\name{permanovas}
\alias{permanovas}
\title{\link{deprecated} Multiple Permanovas on phyloseq (adonis1)}
\usage{
permanovas(
  physeq,
  distm = "bray",
  perms = 999,
  vars = "",
  hellinger = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{physeq}{phyloseq object. A distance matrix is created using \code{otu_table(physeq)}; permanovas are computed for each element of \code{sample_variables(physeq)}}

\item{distm}{choice of dissimilarity index according to \code{\link{vegdist}}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
"gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".}

\item{perms}{the number of permutations to be performed. Default is 999}

\item{vars}{Empty=default tests all metadata variables (should be factor levels >=2). Otherwise give a vector c("var1","var2",...etc)}

\item{hellinger}{boolean. Default TRUE to Hellinger transform the "distance" matrix}

\item{verbose}{boolean. Default is TRUE and prints progress.}
}
\value{
a dataframe containing p-values of each comparison in the first column, and the number of levels for that variable in the second.
}
\description{
Function performs individual, \code{\link{adonis}}-style permanovas for every/selected elements of \code{sample_variables(physeq)}
}
\note{
if \code{sample_data(physeq)} contains NAs, the function will convert these to an additional level "none"

Added the complete adonis output in verbose mode (alex 20200515)

V1 deprecated since we moved to adonis2 having different output (20220610)
}
\examples{

per <- permanovas(phyloseq, "bray")

}
\author{
Alex Bossers \email{a.bossers@uu.nl} inspired by Stephanie Jurburg
}
