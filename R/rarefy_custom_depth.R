#' Convenience function to rarefy a given phyloseq object not to one value for each sample but to allow to rarefy each sample with its own given rarefy target depth.
#'
#' Why? In the exceptional case when you want to proportionally rarefy the data to correct for instance for input volume differences
#' which may be the case in semi-quantitative measurements of the resistome. Since with resistome for example we do not always want to
#' work with relative abundances (100%) and maintaining the ARG total load per sample (FPK or FPKM similar).
#'
#' @param physeq phyloseq object containing at least 2 samples
#' @param rarefy_depth vector of counts per sample to target rarefy (match sample order!)
#' @param seed which random seed to use
#' @return m.final rarefied phyloseq object
#'
#' @author: a.bossers\@uu.nl
#' @note
#'  version: 1.1
#'  versiondate: 20201112
#'
#' @import phyloseq
#' @export


rarefy_custom_depth <- function( physeq, rarefy_depth, seed=2202 ) {
  
  #TODO:
  #FIXME: returns only the last rarefied object within the loop... something with scopes.. :-( Use as single code no function and wors!
  stop("Needs fix and for now run in local code copy mode :-(")
  
  #declare here the function inputs manually as a workaround
#  physeq <- ps                 # phyloseq object.
#  rarefy_depth <- somevector   # vector with value per sample to rarefy to.
#  seed <- 2202                 # just organise the random seed. 
  
  nsamples <- length(sample_names(physeq))

  if( ! nsamples == length(rarefy_depth) ) {
    stop("Exit: Unequal length sample sum vector and N samples in phyloseq object.\nYou have to provide a vector of rarefy_depths having same length as the number of samples.")
  }

  # first do 1st sample then build around it the remaining samples
  # assumes we always have more then 1 sample in the ps object

  i=1
  m.tmp <- subset_samples( physeq, sample_names(physeq) == sample_names(physeq)[i] )
  m.tmp <- rarefy_even_depth( m.tmp, sample.size=rarefy_depth[i], trimOTUs=FALSE, rngseed=seed )
  m.final <- m.tmp

  # now add remainder samples
  for( i in 2:nsamples) {
    m.tmp <- subset_samples( physeq, sample_names(physeq) == sample_names(physeq)[i] )
    m.tmp <- rarefy_even_depth( m.tmp, sample.size=rarefy_depth[i], trimOTUs=FALSE, rngseed=seed )
    m.final <- merge_phyloseq( m.final, m.tmp )
  }

  return( m.final )
}
