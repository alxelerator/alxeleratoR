#' Calculate the co-occurrence of taxa (pair-wise) between one or two tax tables in phyloseq objects
#' 
#' Alex adapted to fit two different matrices to test co-occurrences (i.e. bacteria to resistome)
#'    Also implemented multiple testing corrections
#'    Also allows two differen datatypes over the same samples to be tested
#' 
#' Method details:
#'    Using Spearman's rank (default), pearson or kendall tau as methods.
#' 
#'    Requires one or two phyloseq objects (subset where appropriate first)
#'    	- Providing one ps will calculate WITHIN the ps
#'    	- Providing two PS objects will calculate BETWEEN
#' 
#'    Provide filtering cutoffs for rho and p
#' 
#'    By default it will filter correlations where one or both columns are only zero (filt_zero)
#' 
#'    Verbose just prints a % text line to monitor progress.
#' 
#' @param ps_x phyloseq object X (obligatory)
#' @param ps_y phyloseq object Y (optional)
#' @param method method to use (spearman default). Basically all methods corr.test supports.
#' @param alternative (optional) - alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". "greater" corresponds to positive association, "less" to negative association.
#' @param exact (optional) default FALSE - a logical indicating whether an exact p-value should be computed. Used for Kendall's tau and Spearman's rho.
#' @param adjust default Bejamini-Hochberg (BH). Any method from p.adjust function should work. Provide "none" for raw.
#' @param rho.filt default 0.6
#' @param p.filt filter adjusted p value. default 0.01. if raw p use adjust method "none"
#' @param filt_zero default TRUE. If abundance of taxon in both datasets is zero, that taxon is filtered out and given rho=0 
#' @param verbose default TRUE
#'
#' @returns Co-correlation data frame of all taxa (including self-self)
#'
#' @import phyloseq
#' @export
#'
#' @author Alex Bossers \email{a.bossers@uu.nl} 
#' @author Inspired on Jin Choi's within set algorithm plain R for one set
#' @author Inspired by phylosmith package doing co-occurrence on ps objects but within one object
#' 
ps_co_occurrence <- function( ps_x='', ps_y='', 
                              method="spearman", alternative='two.sided', exact=FALSE, adjust="BH",
                              rho.filt=0.6, p.filt=0.01, filt_zero=TRUE, verbose=TRUE ) 
{
  require(phyloseq)
  
  # check input params
  if ( is(ps_x,"phyloseq") ) {
    a = ps_x
  } else {
    stop("Pyloseq object X required!")
  }
  
  # Check if co-occurrence on itself or on two different tables
  if ( is( ps_y, "phyloseq" ) ) {
    b = ps_y
    n.exclude = 0
    cat("Phyloseq Y provided => all vs all taxa (incl self-self complete matrix)\n")
  } else {
    b = ps_x
    n.exclude = 1 # counter to avoid self-self
    cat("Phyloseq Y NOT provided => self all vs all taxa excluding self-self (half matrix)\n")
  }
  
  a.n <- ntaxa(a)
  b.n <- ntaxa(b)
  
  a.otu <- data.frame( otu_table(a) )
  b.otu <- data.frame( otu_table(b) )
  
  # Check that samples are the same
  if ( ! all( row.names(a.otu) == row.names(b.otu) ) ) {
    stop("Sample rows are not equal")
  }
  
  results <- data.frame()
  ps_y_from <- 1
  verbosecount <- 0

  for( x in 1:(a.n - n.exclude) ) {
    if(n.exclude) {
      ps_y_from <- x + n.exclude
    }
    for( y in ps_y_from:b.n ) {
      
      # track abundance and test if zero
      species1.ab <- sum( a.otu[,x] )
      species2.ab <- sum( b.otu[,y] )
      
      test <- cor.test( a.otu[,x], b.otu[,y], method="spearman", alternative=alternative, exact=exact, na.action=na.rm )
      rho<-test$estimate
      p.value<-test$p.value
      
      if ( filt_zero & ( species1.ab <= 0 | species2.ab <= 0 ) ) {
        rho <- 0
        p.value <- 1
      }
      
      # if results meet filters add to df
      if( abs(rho) >= rho.filt & p.value <= p.filt ) {
        new.row <- data.frame( "taxa.X"=taxa_names(a)[x], "taxa.Y"=taxa_names(b)[y], rho, p.value, species1.ab, species2.ab )
        results <- rbind( results, new.row )
      }
    }
    if(verbose) { 
      verbosecount = verbosecount + 1
      verbperc = verbosecount / a.n * 100
      if( verbperc >= 20 ) {
        verbosecount = 0
        cat( paste( round( x / a.n * 100 ), "% Done of dataset.\n")) 
      }
    }
  }
  
  # apply multiple testing correction
  results$p.adjust <- p.adjust( results$p.value, method=adjust )
  cat("\nMultiple testing correction method '",adjust,"' applied\n")
  #filter again
  results <- results[ results$p.adjust <= p.filt, ]
  
  return(results)
}

