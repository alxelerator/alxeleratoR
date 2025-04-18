% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_co_occurrence.R
\name{ps_co_occurrence}
\alias{ps_co_occurrence}
\title{Calculate the co-occurrence of taxa (pair-wise) between one or two tax tables in phyloseq objects}
\usage{
ps_co_occurrence(
  ps_x = "",
  ps_y = "",
  method = "spearman",
  alternative = "two.sided",
  exact = FALSE,
  adjust = "BH",
  rho.filt = 0.6,
  p.filt = 0.01,
  filt_zero = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{ps_x}{phyloseq object X (obligatory)}

\item{ps_y}{phyloseq object Y (optional)}

\item{method}{method to use (spearman default). Basically all methods corr.test supports.}

\item{alternative}{(optional) - alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". "greater" corresponds to positive association, "less" to negative association.}

\item{exact}{(optional) default FALSE - a logical indicating whether an exact p-value should be computed. Used for Kendall's tau and Spearman's rho.}

\item{adjust}{default Bejamini-Hochberg (BH). Any method from p.adjust function should work. Provide "none" for raw.}

\item{rho.filt}{default 0.6}

\item{p.filt}{filter adjusted p value. default 0.01. if raw p use adjust method "none"}

\item{verbose}{default TRUE}

\item{filt_zeo}{default TRUE}
}
\value{
Co-correlation data frame of all taxa (including self-self)
}
\description{
Alex adapted to fit two different matrices to test co-occurrences (i.e. bacteria to resistome)
Also implemented multiple testing corrections
Also allows two differen datatypes over the same samples to be tested
}
\details{
Method details:
Using Spearman's rank (default), pearson or kendall tau as methods.

Requires one or two phyloseq objects (subset where appropriate first)
- Providing one ps will calculate WITHIN the ps
- Providing two PS objects will calculate BETWEEN

Provide filtering cutoffs for rho and p

By default it will filter correlations where one or both columns are only zero (filt_zero)

Verbose just prints a \% text line to monitor progress.
}
\author{
Alex Bossers \email{a.bossers@uu.nl}

Inspired on Jin Choi's within set algorithm plain R for one set

Inspired by phylosmith package doing co-occurrence on ps objects but within one object
}
