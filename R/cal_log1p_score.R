#' Calculate recalibration matrix for heatmaps
#' Apply transformation function of log avoiding 0 (1p)
#'
#' @param x value
#'
#' @returns tranformed value
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @note
#'   Updates:
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 
cal_log1p_score <- function(x){
  log1p(x)
}
