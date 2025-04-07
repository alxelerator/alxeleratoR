#' Plot calculated procrustes residuals
#'
#' @param df <deprecated> do not use. For backward compatibility only.
#' @param proc_df the procrustus data frame (from vegan)
#' @param main Title to the plot
#' @return Base R plot
#' @note Uses the data from procrustes vegan package
#' @export
#' @author Alex Bossers (UU)

procrustes_error<- function(df="", proc_df, main="Procrustes residuals"){

  res     <- residuals(proc_df)
  qn      <- quantile(res)
  plot(res, type='h', col=cols_t, lwd=2,
       xlab="Index study population",
       ylab="Procrustes residual", main= main)
  abline(h = qn[2:4], lty = c(2, 1, 2))

}
