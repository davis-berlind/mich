#' Summary method for mich.fit objects
#'
#' @description
#' Prints a summary of the resulting fit from calling `mich()`, including the
#' ELBO for the fitted model, estimated change-points with `level`-level
#' credible sets.
#'
#' @param object A mich.fit object. Output of running `mich()` on a numeric vector
#'   or matrix.
#' @param level A scalar. A number between (0,1) indicating the significance
#'   level to construct credible sets at.
#' @param max_length An integer. Detection threshold, see `mich_sets()`. Equal
#'   to `log(T)^2` by default.
#' @param ... Additional arguments to be passed to methods.
#'
#' @method summary mich.fit
#'
#' @return A list. A list of summary quantities including:
#'   * `elbo`: The value of the ELBO for the model.
#'   * `converged`: Indicator for whether the model has converged.
#'   * `level`: The significance level used to construct credible sets.
#'   * `L`,`K`,`J`: The number of mean-scp, var-scp, and meanvar-scp components
#'     included in the model.
#'   * `mean_cp`,`var_cp`,`meanvar_cp`: Lists with `cp`, the estimated
#'     change-points, and `sets`, their corresponding `level`-level credible
#'     sets.
#'
#' @export
#'
summary.mich.fit <- function(object, level = 0.95, max_length = NULL, ...) {
  if (is.matrix(object$y)) {
    T <- nrow(object$y)
    if (is.null(max_length)) max_length <- log(T)^2

    mich_summary <- list()
    mich_summary$mean_cp <- mich_sets(object$pi_bar, max_length, level)
    mich_summary$level <- level
    mich_summary$L <- object$L
    mich_summary$elbo <- max(object$elbo)
    mich_summary$converged <- object$converged

    mcp <- data.frame(change.points = mich_summary$mean_cp$cp,
                      lower = sapply(1:length(mich_summary$mean_cp$cp), function(i) min(mich_summary$mean_cp$sets[[i]])),
                      upper = sapply(1:length(mich_summary$mean_cp$cp), function(i) max(mich_summary$mean_cp$sets[[i]])))
    names(mcp) <- c("change.points", paste0("lower.", level,".credible.set"), paste0("upper.", level,".credible.set"))

    cat("Multivariate MICH Model:\n\n")
    cat(paste0("ELBO: ", mich_summary$elbo, "; Converged: ", mich_summary$converged, "\n\n"))
    cat(paste0("L = ", mich_summary$L, " Mean-SCP Component(s); ", length(mich_summary$mean_cp$cp)," Detected Mean Change-Point(s):\n"))
    print(mcp)
  } else {
    T <- length(object$y)
    if (is.null(max_length)) max_length <- log(T)^2

    mich_summary <- list()
    mich_summary$level <- level
    mich_summary$elbo <- max(object$elbo)
    mich_summary$converged <- object$converged

    cat("Multivariate MICH Model:\n\n")
    cat(paste0("ELBO: ", mich_summary$elbo, "; Converged: ", mich_summary$converged, "\n"))

    if (object$L > 0) {
      mich_summary$L <- object$L
      mich_summary$mean_cp <- mich_sets(object$mean_model$pi_bar, max_length, level)
      cat(paste0("\nL = ", mich_summary$L, " Mean-SCP Component(s); ", length(mich_summary$mean_cp$cp)," Detected Mean Change-Point(s):\n"))

      if (length(mich_summary$mean_cp$cp) > 0) {
        mcp <- data.frame(change.points = mich_summary$mean_cp$cp,
                          lower = sapply(1:length(mich_summary$mean_cp$cp), function(i) min(mich_summary$mean_cp$sets[[i]])),
                          upper = sapply(1:length(mich_summary$mean_cp$cp), function(i) max(mich_summary$mean_cp$sets[[i]])))
        names(mcp) <- c("change.points", paste0("lower.", level,".credible.set"), paste0("upper.", level,".credible.set"))
        print(mcp)
      }
    }
    if (object$K > 0) {
      mich_summary$K <- object$K
      mich_summary$var_cp <- mich_sets(object$var_model$pi_bar, max_length, level)
      cat(paste0("\nK = ", mich_summary$K, " Var-SCP Component(s); ", length(mich_summary$var_cp$cp)," Detected Variance Change-Point(s):\n"))

      if (length(mich_summary$var_cp$cp) > 0) {
        vcp <- data.frame(change.points = mich_summary$var_cp$cp,
                          lower = sapply(1:length(mich_summary$var_cp$cp), function(i) min(mich_summary$var_cp$sets[[i]])),
                          upper = sapply(1:length(mich_summary$var_cp$cp), function(i) max(mich_summary$var_cp$sets[[i]])))
        names(vcp) <- c("change.points", paste0("lower.", level,".credible.set"), paste0("upper.", level,".credible.set"))
        print(vcp)
      }
    }
    if (object$J > 0) {
      mich_summary$J <- object$J
      mich_summary$meanvar_cp <- mich_sets(object$meanvar_model$pi_bar, max_length, level)
      cat(paste0("\nJ = ", mich_summary$J, " MeanVar-SCP Component(s); ", length(mich_summary$meanvar_cp$cp)," Detected Mean-Variance Change-Point(s):\n"))

      if (length(mich_summary$meanvar_cp$cp) > 0) {
        mvcp <- data.frame(change.points = mich_summary$meanvar_cp$cp,
                          lower = sapply(1:length(mich_summary$meanvar_cp$cp), function(i) min(mich_summary$meanvar_cp$sets[[i]])),
                          upper = sapply(1:length(mich_summary$meanvar_cp$cp), function(i) max(mich_summary$meanvar_cp$sets[[i]])))
        names(mvcp) <- c("change.points", paste0("lower.", level,".credible.set"), paste0("upper.", level,".credible.set"))
        print(mvcp)
      }
    }
  }
  invisible(mich_summary)
}
