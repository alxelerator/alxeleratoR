#' Plot the procrustus overlay
#'
#' @param x ordinations
#' @param kind default = 1
#' @param choices Ordination axis (default c(1,2))
#' @param to.target bolean (default TRUE) to select tails=Yrot heads=X
#' @param type points, text or none. Default points (p).
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param main Main plot title.
#' @param len Arrows width. Default 0.05
#' @param cex Text centering. Default 0.7
#' @param arc.col Color of lines. Default "blue".
#' @param ... inherits
#'
#' @return Base R plot
#' @note Uses the data from procrustes vegan package
#' @import vegan
#' @export
#' @author Scavenged from somewhere (you? Contact me for credits)'
#' @author Alex Bossers (UU) \email{a.bossers@uu.nl}


plot_procrustes <- function (x, kind = 1, choices = c(1, 2), to.target = TRUE, type = "p",
                             xlab, ylab, main="Procrustes plot", len = 0.05, cex = 0.7, arc.col="blue",
                             ...)
{
  require(vegan)

  type <- match.arg(type, c("points", "text", "none"))
  if (to.target) {
    tails <- x$Yrot[, choices]
    heads <- x$X[, choices]
  }
  else {
    tails <- x$X[, choices]
    heads <- x$Yrot[, choices]
  }
  if (missing(main))
    main <- "Procrustes plot"
  if (kind <= 1) {
    formals(arrows) <- c(formals(arrows), alist(... = ))
    if (missing(xlab))
      xlab <- paste("Dimension", choices[1])
    if (missing(ylab))
      ylab <- paste("Dimension", choices[2])
    xrange <- range(tails[, 1], heads[, 1]) #These are the ranges of the plot
    yrange <- range(tails[, 2], heads[, 2]) #These are the ranges of the plot
    plot(xrange, yrange, xlab = xlab, ylab = ylab, main = main, arc.col=cols_t,
         type = "n", asp = 1, lwd=5, ...)
    if (kind > 0) {
      abline(v = 0, lty = 3, lwd=1.5) #Abline vertical, changed also ltp
      abline(h = 0, lty = 3, lwd=1.5) #Abline horizontal
      if (ncol(x$rotation) == 2) {
        b <- abs(x$rotation[1, 1])
        EPS <- 1e-08
        if (1 - b > EPS && b > EPS) {
          rot <- x$rotation
          abline(0, rot[1, 2]/rot[1, 1], lty = 1, lwd=2) #Abline cross, changed also ltp
          abline(0, rot[2, 2]/rot[2, 1], lty = 1, lwd=2) #Abline cross
        }
        else {
          abline(v = 0, h = 0, lwd=1.5)
        }
      }
      else {
        Y <- x$Yrot %*% t(x$rotation)
        for (k in 1:ncol(Y)) {
          tmp <- matrix(0, nrow = 2, ncol = ncol(Y))
          tmp[, k] <- range(Y[, k])
          tmp <- tmp %*% x$rotation
          lines(tmp[, choices], lty = 1, lwd=1.5)
          text(tmp[2, choices[1]], tmp[2, choices[2]],
               as.character(k))
        }
      }
      if (type != "none") {
        ow <- options(warn = -1)
        arrows(tails[, 1], tails[, 2], heads[, 1], heads[, 2], lwd=2, col = cols_t, len = len, ...) #Add arrows
        options(ow)
        if (type == "text" && !is.null(rownames(tails)))
          ordilabel(tails, cex = cex, pch=19, col = cols_t, ...)
        else points(tails, cex = cex, pch=19, col = cols_t, ...)
      }
    }
    out <- list(heads = heads, points = tails)
    class(out) <- "ordiplot"
  }
  else if (kind == 2) {
    if (missing(xlab))
      xlab <- "Index"
    if (missing(ylab))
      ylab <- "Procrustes residual"
    res <- residuals(x)
    q <- quantile(res)
    plot(res, type = "h", xlab = xlab, ylab = ylab, main = main,
         ...)
    abline(h = q[2:4], lty = c(2, 1, 2))
    out <- list(sites = cbind(seq(along = res), res))
    class(out) <- "ordiplot"
  }
  invisible(out)
}
