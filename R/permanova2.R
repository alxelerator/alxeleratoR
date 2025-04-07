#' Single (multivariable) PERMANOVA v2
#' Uses vegan::adonis2.
#' This function is just a convenience wrapper for vegan::adonis2. Please give proper credits to those authors.
#'
#' Function performs individual, \code{\link{adonis2}} permanova for given variable or (multivariable) formula
#'
#' @param physeq phyloseq object. A distance matrix is created using \code{otu_table(physeq)}; permanovas are computed for each element of \code{sample_variables(physeq)}
#' @param formula Defaults to distnmatrix ~ var
#' @param var ONE variable. Required if no formula is specified! (should be factor levels >=2).
#' @param distm choice of dissimilarity index according to \code{\link{vegdist}}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param perms the number of permutations to be performed. Default is 999
#' @param hellinger TRUE/FALSE. By default performs a Hellinger transformation of the OTU table. Here you can turn it off.
#' @param rngseed set.seed number to use for reproducible permutations. Default rngseed=2202.
#' @param ... Allows additional arguments to the adonis2 function to be taken along. For instance parallel = parallel::detectCores() - 1
#'
#' @return vegan2 default output.
#'
#' @export
#' @import vegan
#' @import phyloseq
#' @import stringr
#' @import dplyr
#'
#' @note if \code{sample_data(physeq)} contains NAs, the function will convert these to an additional level "none"
#'
#' @author Convenience wrapper: Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @examples
#'   per <- permanova2( physeq, distm="bray", formula="~ Treatment + AnimalNumber", var="", rngseed="2202", parallel=(parallel::detectCores() - 1) )
#'   per <- permanova2( physeq, distm="bray", formula="", var="Treatment", rngseed="2202")
#'
permanova2 <- function(physeq="", distm="bray", formula="", var="", permutations=999, hellinger=TRUE, rngseed=2202, ...) {

  require(phyloseq)
  require(vegan)
  require(stringr)
  require(dplyr)
  
  if( class(physeq)[1] != "phyloseq" ) {
    stop("Required phyloseq object is missing.")
  }
  if( var[1] == "" & formula == "" ) {
    stop('Neither a formula or variable to use is given. Specify at least one of these two.\n  Example: permanova2( ps, distm="bray", formula="", var="Treatment")\n  Formula only required for multivariable analysis (defaults to distn ~ VAR).')
  } else if( var[1] != "" & formula != "" ) {
    stop('Both a formula AND variable are given. Specify only one of these two.\n  Example: permanova2( ps, distm="bray", formula="", var="Treatment")\n  Formula only required for multivariable analysis (defaults to distn ~ VAR).')
  }

  if( formula != "" ) {
    #remove any spaces
    formula <- gsub( " ", "", formula, fixed=TRUE )
    #add ~ if not present
    if(substr(formula,1,1)!="~") {
      formula <- paste0("~",formula)
    }
    message(paste("Final formula:", formula))
    #check if variables exist and throw warning since not sure we capture each possibility of splitting here
    vars <- str_split(formula, "[~+*]+")[[1]] #was "[ ~ +*|()0-9]+")[[1]]
    vars <- vars[ vars > "" ]
    if( length(vars) != length( intersect( vars, sample_variables(physeq) ) ) ){
      warning(paste("Variables: ",length(vars),"\n"))
      warning("Seems like one or more of the specified variables in your formula are not present in the metadata? \n    - Check spelling and UpperLowerCase.\n    - Does formula start with ~?")
    }
  }

  if( length( var[ var > "" ] ) >= 1 ) {
    #check if variable exist and is only length 1
    if( length(var) > 1 ){
      stop("More then one variable is specified!\n  Use permanovas2() {note the s} for testing multiple variables at once!\n  Or specify a formula containing these vars in multivariable setting.")
    } else if( ! var %in% sample_variables(physeq) ) {
      stop("Specified variable is not present in the metadata. Check spelling and UpperLowerCase.")
    }
    vars <- var
  }

  if (physeq@otu_table@taxa_are_rows == TRUE ) {
    otutab <- t( otu_table(physeq) )
  } else {
    otutab <- ( otu_table(physeq) )
  }
  
  if(hellinger) {
    message( paste0("Using Hellinger transformed data and distance metric '",distm,"'") )
    distn <- vegdist(decostand(otutab, "hellinger"), method=distm)
  } else {
    message( paste0("Using plain (not Hellinger transformed) data and distance metric '",distm,"'") )
    distn <- vegdist( otutab, method=distm )
  }

  if( formula != "" ){
    formula = paste("distn",formula)
  } else {
    formula = paste("distn ~",var)
  }
  formula_to_use <- as.formula(formula)
  
  meta.data <- data.frame( sample_data(physeq) ) %>%
    select( vars )

  #check metadata if NAs and levels > 1
  for (i in seq_along (vars)) {
    if ( anyNA( meta.data[[i]] ) ){
      levels( meta.data[[i]] ) <- c( levels( meta.data[[i]] ), "none")
      meta.data[[i]][ is.na( meta.data[[i]] ) ] <-  "none"
      warning( paste0("<NA>s in metadata column '",vars[i],"' have been converted to an additional category 'none'") )
    }
    if ( nlevels( meta.data[[i]] ) == 1 ) {
      stop( paste0( "Variable '",vars[i],"', has only one level!") )
    }
  }
  
  #calc
  set.seed(2202)
  adonisResult <- adonis2( formula_to_use, data=meta.data, ... )

  return(adonisResult)
}


