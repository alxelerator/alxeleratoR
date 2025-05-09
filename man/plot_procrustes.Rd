% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_procrustes.R
\name{plot_procrustes}
\alias{plot_procrustes}
\title{Plot calculated procrustus overlay}
\usage{
plot_procrustes(
  x,
  kind = 1,
  choices = c(1, 2),
  to.target = TRUE,
  type = "p",
  xlab,
  ylab,
  main = "Procrustes plot",
  len = 0.05,
  cex = 0.7,
  arc.col = "blue",
  ...
)
}
\arguments{
\item{x}{ordinations}

\item{kind}{default = 1}

\item{choices}{Ordination axis (default c(1,2))}

\item{to.target}{bolean (default TRUE) to select tails=Yrot heads=X}

\item{type}{points, text or none. Default points (p).}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{main}{Main plot title.}

\item{len}{Arrows width. Default 0.05}

\item{cex}{Text centering. Default 0.7}

\item{arc.col}{Color of lines. Default "blue".}

\item{...}{inherits}
}
\value{
Base R plot
}
\description{
Plot calculated procrustus overlay
}
\note{
Uses the data from procrustes vegan package
}
\author{
Scavenged from somewhere (you? Contact me for credits)'

Alex Bossers (UU) \email{a.bossers@uu.nl}
}
