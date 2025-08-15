#' Plot method for mich.fit objects
#'
#' @description
#' Plots the resulting fit and estimated change-points from calling `mich()`.
#'
#' @param x A mich.fit object. Output of running `mich()` on a numeric vector
#'   or matrix.
#' @param level A scalar. A number between (0,1) indicating the significance
#'   level to construct credible sets at when `cs == TRUE`.
#' @param max_length An integer. Detection threshold, see `mich_sets()`. Equal
#'   to `log(T)^2` by default.
#' @param signal A logical. If `TRUE`, then the posterior mean and precision
#'   signals are also plotted.
#' @param cs A logical. If `TRUE`, then `level`-level credible sets for each
#'   detected change-point are also plotted.
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
#' @examples
#' set.seed(222)
#' # generate univariate data with two mean-variance change-points
#' y = c(rnorm(100,0,10), rnorm(100,10,3), rnorm(100,0,6))
#' fit = mich(y, J = 2) # fit two mean-variance change-points
#' # plot change-points with 95% credible sets
#' plot(fit, level = 0.95, cs = TRUE)
#' # fit one mean and one mean-variance change-point
#' fit = mich(y, J = 1, L = 1)
#' # plot change-points with 95% credible sets and signal
#' plot(fit, level = 0.95, cs = TRUE, signal = TRUE)
#'
#' # generate correlated mulitvariate data with two mean-variance change-points
#' T <- 150
#' Sigma <- rbind(c(1, 0.7), c(0.7, 2))
#' d <- ncol(Sigma)
#' Sigma_eigen <- eigen(Sigma)
#' e_vectors <- Sigma_eigen$vectors
#' e_values <- Sigma_eigen$values
#' Sigma_sd <- e_vectors %*% diag(sqrt(e_values)) %*% t(e_vectors)
#' Z <- sapply(1:d, function(i) rnorm(T))
#' mu <- c(-1, 2)
#' mu_t <- matrix(0, nrow = 70, ncol=d)
#' mu_t <- rbind(mu_t, t(sapply(1:30, function(i) mu)))
#' mu_t <- rbind(mu_t, matrix(0, nrow = 50, ncol = d))
#' Y <- mu_t + Z %*% Sigma_sd
#' # fit MICH and pick L automatically using ELBO
#' fit = mich(Y, L_auto = TRUE)
#' plot(fit, level = 0.95, cs = TRUE, signal = TRUE)
#'
plot.mich.fit <- function(x, level = 0.95, max_length = NULL, signal = FALSE, cs = TRUE, n_plots = 5, ...) {
  fit <- x
  # multivariate plots
  if (is.matrix(fit$y)) {
    par(ask = TRUE)
    T <- nrow(fit$y)
    if (is.null(max_length)) max_length <- log(T)^2

    if (fit$L > 0) {
      N <- ncol(fit$y) %/% n_plots
      M <- ncol(fit$y) %% n_plots
      cp_est <- mich_sets(fit$pi_bar, max_length, level)

      for (i in seq_len(N)) {
        par(mfrow = c(n_plots, 1),
            mar = c(4,4,4,3),
            oma = c(0,0,0,0))
        for (j in 1:n_plots) {
          if (j > 1) {
            par(mar = c(4,4,0.5,3))
          }
          plot(fit$y[, (i-1) * n_plots + j],
               main = ifelse(j==1, "Mean Change-Points",""),
               xlab = ifelse(j == n_plots, "Index", ""),
               ylab = ifelse(is.null(names(fit$y[, (i-1) * n_plots + j])),
                             paste0("Series ",  (i-1) * n_plots + j),
                             names(fit$y[, (i-1) * n_plots + j])),
               type = "l")
          if (signal) {
            lines(fit$mu[, (i-1) * n_plots + j], col = "blue", lwd = 2)
          }
          if (cs) {
            for (k in unlist(cp_est$sets)) {
              rect(xleft = k - 0.5, xright = k+ 0.5,
                   ybottom = par("usr")[3], ytop = par("usr")[4],
                   col =  adjustcolor("lightblue", alpha.f = .8),
                   border = NA)
            }
          }
          abline(v = cp_est$cp, lty = 2, col = "red", lwd = 2)
        }
      }
      if (M > 0) {
        par(mfrow = c(M, 1),
            mar = c(4,4,4,3),
            oma = c(0,0,0,0))
        for (j in 1:M) {
          if (j > 1) {
            par(mar = c(4,4,0.5,3))
          }
          plot(fit$y[, N * n_plots + j],
               main = ifelse(j==1, "Mean Change-Points",""),
               xlab = ifelse(j == M, "Index", ""),
               ylab = ifelse(is.null(names(fit$y[, N * n_plots + j])),
                             paste0("Series ",  N * n_plots + j),
                             names(fit$y[, N * n_plots + j])),
               type = "l")
          if (signal) {
            lines(fit$mu[, N * n_plots + j], col = "blue", lwd = 2)
          }
          if (cs) {
            for (k in unlist(cp_est$sets)) {
              rect(xleft = k - 0.5, xright = k + 0.5,
                   ybottom = par("usr")[3], ytop = par("usr")[4],
                   col =  adjustcolor("lightblue", alpha.f = .8),
                   border = NA)
            }
          }
          abline(v = cp_est$cp, lty = 2, col = "red", lwd = 2)
        }
      }
    }
    par(ask = FALSE)
  } else {
    # univariate plots
    T <- length(fit$y)
    if (is.null(max_length)) max_length <- log(T)^1.5
    par(mfrow = c((fit$L > 0) + (fit$K > 0) + (fit$J > 0), 1),
        mar = c(4,4,4,3),
        oma = c(0,0,0,0))
    if (fit$L > 0) {
      cp_est <- mich_sets(fit$mean_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Mean Change-Points",  ylab = "y")
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 2)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
      }
      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i - 0.5, xright = i + 0.5,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 2)
    }
    if (fit$K > 0) {
      cp_est <- mich_sets(fit$var_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Variance Change-Points",  ylab = "y")
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 2)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
      }
      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i-0.5, xright = i + 0.5,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 2)
    }
    if (fit$J > 0) {
      cp_est <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
      plot(fit$y, type = "l", main = "Mean-Variance Change-Points", ylab = "y")
      if (signal) {
        lines(fit$mu, col = "blue", lwd = 2)
        lines(fit$mu+2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
        lines(fit$mu-2*1/sqrt(fit$lambda), col = "blue", lwd = 2, lty = 2)
      }
      if (cs) {
        for (i in unlist(cp_est$sets)) {
          rect(xleft = i-0.5, xright = i + 0.5,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col =  adjustcolor("lightblue", alpha.f = .8),
               border = NA)
        }
      }
      abline(v = cp_est$cp, lty = 2, col = "red", lwd = 2)
    }
  }
}
