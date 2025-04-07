#' Initialise environment
#' Personal preferences to set for instance cowplot theme globally for ggplot2
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @import cowplot
#' @export
#' @note
#'   Updates:
#'          2020-02-03 Alx initial checkin
#' 
alxinit <- function(){
  
  # theme
  suppressWarnings(library(cowplot))
  theme_set(theme_cowplot())

  # convenience
  today <- lubridate::today()

  # find home
  if(.Platform$OS.type == "unix") {
    home <- "~"
  } else {
    home <- "D:"
  }

  #seed
  set.seed(2202)
}
