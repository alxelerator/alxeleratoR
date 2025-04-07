#' Apply the centered log-ratio trasformation on a relative abundance matrix!
#' Species in columns and samples in rows.
#' Defined as log^base( x / geometric_mean(rowx) )
#'
#' @param x OTU/ASV relative abundance matrix species as cols and samples as rows.
#' @param base default=base log 2 (can be 10 or exp(1)....
#' @param pseudocount default=FALSE applies a pseudocount of MIN(relative abundance)/2 to exact zero relative abundance entries in ASV/OTU table before taking logs.
#' @return transformed matrix
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @author joey711 Original idea scavenged from R phyloseq Shiny documentation and slightly adapted https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
#'
#' Versions:
#'    20211103 mad pseudocount non default (was not intuitive)
#'    20200909 added default pseudocount to prevent log(0)
#'    20190322 initial scavenge and adaptation'
#' @export

clr_transform <- function( x, base=2, pseudocount=FALSE )
{
	
	gm_mean = function(x, na.rm=TRUE)	{
  		# The geometric mean, with some error-protection bits.
  		exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
	}

	# First transform zero's using pseudocount MIN(rel_abundance)/2
	if( pseudocount ) {
	  # determine the minimum rel abundance in the table (larger then 0)
	  pc <- min( x[ x>0 ] ) / 2
	  x[ x==0 ] <- pc
	}

  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0

  return(x)
}
