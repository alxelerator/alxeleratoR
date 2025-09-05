#' Calculate the prevalence at lowest taxonomic rank (i.e. ASV) from phyloseq. Function is used in the prevalence_plot() function
#'
#' Function will return a prevalence data frame from each OTU
#' N occurences in samples (y) and total abundance (x)
#' Optional the prevalence can be relative (abundance divided by N samples)'
#'
#' @param ps a phyloseq object containing cummunity abundance data
#' @param prev_relative prevalence as numbers of samples or relative to total number of samples
#' @return dataframe object
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @author Alex Bossers \email{alex.bossers@wur.nl}
#''
#' @note Prevalence plot based on Callahan et al 2016 (McMurdie and Holmes group) https://f1000research.com/articles/5-1492/v2
#'
#' @export
#'
#' @examples
#' # prev <- prevalence_taxa( ps = PhyObject )
#'
#' @import phyloseq

prevalence_taxa <- function( ps, prev_relative = FALSE ) {
    require(phyloseq)
    # variable rank deprecated.
    if ( prev_relative ) {
        correctionfactor = nsamples( ps )
    } else {
        correctionfactor = 1
    }
    prev0  <- apply(X = otu_table( ps ),
                    MARGIN = ifelse(taxa_are_rows( ps ), yes = 1, no = 2),
                    FUN = function( x ) { sum( x > 0 ) / correctionfactor } )
    prevdf <- data.frame( Prevalence = prev0, TotalAbundance = taxa_sums( ps ), tax_table( ps ))
    return(prevdf)
}

