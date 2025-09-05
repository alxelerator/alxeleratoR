#' Multiple Permanovas v2
#' Uses vegan::adonis2 which returns a slight different object.
#' V2 now returns a proper data frame of tested vars and R2 and F.
#'
#' Function performs individual, adonis-style permanovas for every/selected elements of sample_variables(physeq)
#'
#' @param physeq phyloseq object. 
#' @param distm choice of dissimilarity index according to vegan::vegdist. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
#'              "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param perms the number of permutations to be performed. Default is 999
#' @param vars Empty=default tests all metadata variables (should be factor levels >=2). Otherwise give it a vector c("var1","var2",...etc) of factor variables.
#' @param rngseed set.seed number to use for reproducible permutations. Default rngseed=2202
#' @param verbose boolean. Default is TRUE and prints progress.
#' @param hellinger boolean. Default TRUE to Hellinger transform the "distance" matrix
#'
#' @return a dataframe containing p-values of each comparison in the first column, and the number of levels for that variable in the second.
#'
#' @export
#' @import vegan
#' @import phyloseq
#' @note if sample_data(physeq) contains NAs, the function will convert these to an additional level "none"
#' @note Added the complete adonis output in verbose mode (alex 20200515)
#' @note V2 Moved to adonis2 since adonis will be deprecated soon (20220610)
#'
#' @author Alex Bossers \email{a.bossers@uu.nl} based on initial idea by Stephanie Jurburg.
#' @examples
#' # per <- permanovas2(ps, distm="bray")
#'
permanovas2 <- function(physeq, distm="bray", perms=999, vars="", hellinger=TRUE, rngseed=2202, verbose=TRUE ) {

  require(phyloseq)
  require(vegan)

  if( vars[1] == "" ) {
    vars <- as.data.frame(lapply(as.data.frame(as.matrix(sample_data(physeq))), factor))
  } else {
    vars <- as.data.frame(lapply(as.data.frame(as.matrix(sample_data(physeq)[,vars])), factor))
  }

  ADONIS <- matrix(nrow=length(vars), ncol=4)
  row.names(ADONIS) <- colnames(vars)
  colnames(ADONIS) <- c( "Levels", "R2", "F", "Pr(>F)" )
  if (physeq@otu_table@taxa_are_rows == TRUE ) {
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

  #loop over vars
  for (i in seq_along (vars)) {
    if (anyNA(vars[[i]])){
      levels(vars[[i]]) <- c(levels(vars[[i]]), "none")
      vars[[i]][is.na(vars[[i]])] <-  "none"
    }
    if ( nlevels(vars[[i]]) == 1 ) {
      message( paste0( "Skipping variable '",rownames(ADONIS)[i],"', not enough levels!") )
      next
    }

    #make reproducible
    set.seed(rngseed)
    
    adonisResult <- adonis2(distn~as.factor(na.pass(vars[[i]])), permutations=perms)
    if(verbose){
      cat(paste0("\nVariable: ",colnames(vars)[i]),"\n")
      print(adonisResult)
      cat("\n\n")
    }
    ADONIS[i,"Levels"] <- nlevels( vars[[i]] )
    ADONIS[i,"R2"] <- adonisResult[1,"R2"]
    ADONIS[i,"F"] <- adonisResult[1,"F"]
    ADONIS[i,"Pr(>F)"] <- adonisResult[1,"Pr(>F)"]
  }
  return(ADONIS)
}


