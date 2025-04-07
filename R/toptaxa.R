#' Returns the top-x list of taxa names at given rank for given samples based on its relative abundances.
#'
#' @param ps phyloseqobject (does not have to be relative abundance, will be converted on the fly)
#' @param rank optinal taxonomic rank (default="Phylum").
#' @param top optional top number of taxa to return (default 10)
#' @param samples optional name of sample or set of sample names to restrict to (default NA, all samples)
#' @return toplist taxa
#' @import phyloseq
#'
#' @export
#' @author Alex.Bossers@wur.nl
#' Version 1: 20190424'


toptaxa <- function( ps, rank="Phylum", top=10, samples=NA, NArm=FALSE ) {
  require(phyloseq)
  if( ! is.na(samples[1]) ) {
    pstemp = prune_samples( samples, ps)
  } else {
    pstemp = ps
  }
  pstemp = tax_glom2( pstemp, taxrank = rank, NArm=NArm )
  pstemp = transform_sample_counts( pstemp, function(x) 100*x/sum(x) )

  if( ntaxa(pstemp) < top ) {
    top = ntaxa(pstemp)
    warning("Number of taxa in object is less then requested top list! Reduced to ", top)
  }

  toplist = names( sort( taxa_sums(pstemp), decreasing = TRUE ) )[1:top]
  return( toplist )
}
