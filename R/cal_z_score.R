#' Calculate recalibration matrix for heatmaps
#' Apply transformation function
#' Z score is difference to mean divided by sd
#'
#' @param x value
#'
#' @returns tranformed value
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @export
#' 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
