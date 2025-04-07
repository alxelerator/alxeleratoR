#' Calculate recalibration matrix for heatmaps
#' Apply transformation function
#' Z score is difference to mean divided by sd
#'
#' @param x value
#'
#' @returns tranformed value
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @note
#'   Updates:
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
