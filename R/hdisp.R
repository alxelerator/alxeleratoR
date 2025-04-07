#' Multiple tests for homogeneity of dispersions (vegan::betadisper)
#'
#' Function performs individual, \code{\link{betadisper}}-type tests for homogeneity of dispersions for every element of \code{sample_variables(physeq)}. Community matrix is Hellinger-transformed prior to the calculation of distances using specified matrix.
#'
#' @param physeq a \link{phyloseq} object
#' @param distm choice of dissimilarity index according to \link{vegdist}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param vars selected variables from the ps object (should be type factor and 2 or more levels). Default uses all.
#' @param hellinger default=TRUE. Performs Hellinger transformation on the otu table.
#' 
#' @return function prints a list of p-values and returns a list of \code{betadisper} outputs which can be examined with \code{output[[x]]} to view the output for the xth element of \code{sample_data(physeq)}
#'
#' @import phyloseq
#' @import vegan
#' @import foreach
#'
#' @export
#'
#' @examples
#' q <- hdisp(physeq, distm="bray", vars=c("Type","Date","Health"))
#' plot(q[[1]])
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#'
#' Alx: added return list incl p-values
#'
hdisp <- function(physeq, distm="bray", vars=NULL, hellinger=TRUE) {
  require(foreach)
  require(phyloseq)
  require(vegan)
  if( length(vars) > 0 ) {
    #subset the variables from metadata
    cat("Subsetting selected metadata variables!\n")
    physeq@sam_data <- physeq@sam_data[,vars]
  } else {
    cat("Using all metadata variables!\n")
  }

  if (physeq@otu_table@taxa_are_rows == TRUE ) {
    otutab <- t( otu_table(physeq) )
  } else {
    otutab <- ( otu_table(physeq) )
  }
  
  if(hellinger) {
    message( paste0("Using Hellinger transformed data and distance metric '",distm,"'") )
    dist.mat <- vegdist(decostand(otutab, "hellinger"), method=distm)
  } else {
    message( paste0("Using plain (not Hellinger transformed) data and distance metric '",distm,"'") )
    dist.mat <- vegdist( otutab, method=distm )
  }
  
  vars <- as.data.frame(lapply(as.data.frame(as.matrix(sample_data(physeq))), factor))
  bdisp <- matrix(nrow=length(vars), ncol=1)
  row.names(bdisp) <- colnames(vars)
  betadisper.results <- foreach(i=1:length(vars)) %do% betadisper(dist.mat,vars[,i])

  for (i in 1:length(vars)) {
    bdisp[i] <- anova(betadisper.results[[i]])$`Pr(>F)`[1]
  }
  colnames(bdisp)="p.values"
  #print(bdisp)
  return( list( Pvalues = bdisp, betadisper = betadisper.results ) )
}

