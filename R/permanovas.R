#' [deprecated] Multiple Permanovas on phyloseq (adonis1)
#'
#' Function performs individual, \code{\link{adonis}}-style permanovas for every/selected elements of \code{sample_variables(physeq)}
#'
#' @param physeq phyloseq object. A distance matrix is created using \code{otu_table(physeq)}; permanovas are computed for each element of \code{sample_variables(physeq)}
#' @param distm choice of dissimilarity index according to \code{\link{vegdist}}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
#'              "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param perms the number of permutations to be performed. Default is 999
#' @param vars Empty=default tests all metadata variables (should be factor levels >=2). Otherwise give a vector c("var1","var2",...etc)
#' @param hellinger boolean. Default TRUE to Hellinger transform the "distance" matrix
#' @param verbose boolean. Default is TRUE and prints progress.
#'
#' @return a dataframe containing p-values of each comparison in the first column, and the number of levels for that variable in the second.
#'
#' @export
#' @import vegan
#' @import phyloseq
#' @note if \code{sample_data(physeq)} contains NAs, the function will convert these to an additional level "none"
#' @note Added the complete adonis output in verbose mode (alex 20200515)
#' @note V1 deprecated since we moved to adonis2 having different output (20220610)
#'
#' @author Alex Bossers \email{a.bossers@uu.nl} inspired by Stephanie Jurburg
#' @examples
#'
#' per <- permanovas(phyloseq, "bray")
#'
permanovas <- function(physeq, distm="bray", perms=999, vars="", hellinger=TRUE, verbose=TRUE) {

  stop("Please migrate to permanovas2 since vegan::adonis will be deprecated and replaced with adonis2 that has slightly different output!")
  
  if( vars[1] == "" ) {
    vars <- as.data.frame(lapply(as.data.frame(as.matrix(sample_data(physeq))), factor))
  } else {
    vars <- as.data.frame(lapply(as.data.frame(as.matrix(sample_data(physeq)[,vars])), factor))
  }

  ADONIS <- matrix(nrow=length(vars), ncol=2)
  row.names(ADONIS) <- colnames(vars)
  colnames(ADONIS) <- c("P.value", "Levels")
  if (physeq@otu_table@taxa_are_rows ==TRUE ) {
    otutab <- t(otu_table(physeq))
  } else {
    otutab <- (otu_table(physeq))
  }

  if(hellinger) {
    message( paste0("Using Hellinger transformed data and distance metric '",distm,"'") )
    distn <- vegdist(decostand(otutab, "hellinger"), method=distm)
  } else {
    message( paste0("Using plain (not Hellinger transformed) data and distance metric '",distm,"'") )
    distn <- vegdist( otutab, method=distm )
  }

  for (i in seq_along (vars)) {
    if (anyNA(vars[[i]])){
      levels(vars[[i]]) <- c(levels(vars[[i]]), "none")
      vars[[i]][is.na(vars[[i]])] <-  "none"
    }
    if ( nlevels(vars[[i]]) == 1 ) {
      message( paste0( "Skipping variable '",rownames(ADONIS)[i],"', not enough levels!") )
      next
    }

    adonisResult <- adonis(distn~as.factor(na.pass(vars[[i]])), permutations=perms)
    if(verbose){
      cat(paste0("\nVariable: ",colnames(vars)[i]),"\n")
      print(adonisResult)
      cat("\n\n")
    }
    ADONIS[i,1] <- adonisResult$aov.tab$`Pr(>F)`[1]
    ADONIS[i,2] <- nlevels(vars[[i]])
  }
  return(ADONIS)
}


