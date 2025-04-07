#' Generates a distinct scale of n=x elements
#' based on QUAL Rcolorebrewer set
#' Manually replaced a few too dominant or too dark colors with other variants
#' 
#' @param n number of desired distinct colors
#' @param lightfirst boolean. Default TRUE.
#'
#' @returns vector of hex colors size N
#' @import RColorBrewer
#' @import pals
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @note Will modify this function to have repeated HUE variants of other frequently used colorschemes
#'
#' @note
#'   Updates:
#'          2021-12-02 Alx update two subsequent grey colors to more blueish
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 
colorDistinct <- function(n,lightfirst=TRUE) {

  require(RColorBrewer)
  require(pals)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #change  a few too bright colors
  col_vector[6] <- "#F60A7F"
  col_vector[22] <- "#990000" #brightred to darkred
  col_vector[28] <- "#A14918"
  col_vector[37] <- "#d0d0FF" # prevent a subsequent double light grey 
  if( lightfirst ) { col_vector <- rev(col_vector) }  #lighter colors first?

  #enlarge palette by repeating 4 times the same palet. Not useful for color bu prevents errors on large sets.
  col_vector <- c( col_vector, col_vector, col_vector, col_vector )

  return( head(col_vector,n) )
}
