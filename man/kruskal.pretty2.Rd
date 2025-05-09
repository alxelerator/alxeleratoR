% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kruskal_pretty2.R
\name{kruskal.pretty2}
\alias{kruskal.pretty2}
\title{Calculates the Kruskal-Wallis rank-sum test on SIMPER pretty output}
\usage{
kruskal.pretty2(ps, csv, interesting, output_name, taxrank = "")
}
\arguments{
\item{ps}{phyloseq object}

\item{csv}{CSV output file from simper.pretty2}

\item{interesting}{columns of var of interest in metadata (same as simper.pretty input)}

\item{output_name}{desired name of output (i.e. outputname_krusk_simper.csv)}

\item{taxrank}{For example Kingdom or Genus. If nothing/empty string given, a selected Taxonomic concat will be made}
}
\value{
CSV df (and file saved)
}
\description{
Debugged, changed and extended to allow an phyloseq object and tax rank selection
}
\details{
What it does:
This script takes the output .csv of simper_pretty.R, and the OTU/metadata/taxonomy tables, and performs the non-parametric Kruskal-Wallis rank-sum test
on each OTU in the .csv file. Output is a .csv file containing the same contents of simper.pretty output with the following info: p-value, fdr corrected
p-value, OTU taxonomic classification (if applicable), mean rel. abund and std dev of otu/tax_lvl in group 1 of comparison, and mean rel. abund and std
dev of otu/tax_lvl in group 2 of comparison.
}
\note{
csv = clean_simper.csv (*Must be imported as data.frame; i.e. csv= data.frame(read.csv("PATH to .csv")))

kruskal.pretty: automates kruskal test execution on simper_pretty.R output
}
\author{
Andrew Steinberger asteinberger@wisc.edu   Suen Lab, University of Wisconsin-Madison

Alex Bossers \email{a.bossers@uu.nl} for phyloseq and extensions

\if{html}{\out{<div class="sourceCode">}}\preformatted{  Copyright (C) 2016 Andrew Steinberger

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
}\if{html}{\out{</div>}}

Changelog

###########################################################################################
Example of use
kruskal.pretty(bacteria, metrics, csv, c('sloth_sp','type','sp.type'), 'sloth', taxonomy)
}
