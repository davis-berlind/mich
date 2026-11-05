#' Plot method for mich.fit objects
#'
#' @description
#' Plots the resulting fit and estimated change-points from calling `mich()`.
#'
#' @param x A mich.fit object. Output of running `mich()` on a numeric vector
#'   or matrix.
#' @param level A scalar. A number between (0,1) indicating the significance
#'   level to construct credible sets at when `cs == TRUE`.
#' @param signal A logical. If `TRUE`, then the posterior mean and precision
#'   signals are also plotted.
#' @param cs A logical. If `TRUE`, then `level`-level credible sets for each
#'   detected change-point are also plotted.
#' @param detect A logical. If `TRUE`, then credible sets are only constructed
#'   for detected changes.
#' @param n_plots An integer. Number of plot to display at once when data `y` is
#'   a matrix.
#' @param ... Additional arguments to be passed to methods.
#'
#' @method plot mich.fit
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline lines par rect
#'
#' @return Invisibly returns `NULL`.
#' @export
#'
plot.mich.fit <- function(
    x,
    level = 0.95,
    signal = FALSE,
    cs = TRUE,
    detect = TRUE,
    n_plots = 5,
    ...
  ) {
  fit <- x
  # multivariate plots
  if (is.matrix(fit$y)) {
    T <- nrow(fit$y)
    d <- ncol(fit$y)
    n_plots <- min(d, n_plots)

    if (fit$L > 0) {
      N <- ceiling(d / n_plots)
      if (N > n_plots) par(ask = TRUE)
      cp_est <- mich_sets(fit$pi_bar, level, detect)

      for (i in seq_len(N)) {
        n_plots <- ifelse(i <= d %/% n_plots, n_plots, d %% n_plots)
        par(mfrow = c(n_plots, 1),
            #mar = c(4,4,4,3),
            oma = c(0,0,0,0))
        for (j in 1:n_plots) {
          if (j > 1) {
            #par(mar = c(4,4,0.5,3))
          }
          plot(fit$y[, (i-1) * n_plots + j],
               main = ifelse(j==1, "Mean Change-Points",""),
               xlab = ifelse(j == n_plots, "Index", ""),
               ylab = ifelse(is.null(names(fit$y[, (i-1) * n_plots + j])),
                             paste0("Series ",  (i-1) * n_plots + j),
                             names(fit$y[, (i-1) * n_plots + j])),
               type = "l", ... = ...)

          if (cs) {
            for (k in unlist(cp_est$sets)) {
              rect(xleft = k - 1, xright = k + 1,
                   ybottom = par("usr")[3], ytop = par("usr")[4],
                   col =  adjustcolor("lightblue", alpha.f = .8),
                   border = NA)
            }
          }
          if (signal) {
            lines(fit$mu_bar[, (i-1) * n_plots + j], col = "blue", lwd = 1)
          }
          abline(v = cp_est$cp, lty = 2, col = "red", lwd = 1)
        }
      }
    }
    par(ask = FALSE)
  } else {
    # univariate plots
    T <- length(fit$y)
    if (is.null(max_length)) max_length <- log(T)^1.5
    par(mfrow = c((fit$L > 0) + (fit$K > 0) + (fit$J > 0), 1))
    if (fit$L > 0) {
      cp_est <- mich_sets(fit$mean_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Mean Change-Points",  ylab = "y", ... = ...)

      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i - 1, xright = i + 1,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 1)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 1)
    }
    if (fit$K > 0) {
      cp_est <- mich_sets(fit$var_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Variance Change-Points",  ylab = "y", ... = ...)

      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i - 1, xright = i + 1,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 1)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 1)
    }
    if (fit$J > 0) {
      cp_est <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Mean-Variance Change-Points", ylab = "y", ... = ...)

      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i - 1, xright = i + 1,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 1)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 1, lty = 2)
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 1)
    }
  }
}
