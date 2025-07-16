#' Log Weighted Mean-SCP Prior
#'
#' @description
#' Log weighted prior for the Mean-SCP model as described in Appendix C.2 of
#' Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that in the
#' absence of a change-point, the posterior probabilities are approximately
#' uniform.
#'
#'
#' @param T An integer. Number of observations (rows) of \eqn{\mathbf{y}_{1:T}}.
#' @param d An integer. Dimension (columns) of \eqn{\mathbf{y}_{1:T}}.
#'
#' @return A numeric vector. Log prior probabilities.
#' @export
#'
log_mean_prior <- function(T, d) {
  log_pi <- rep(0, T)
  for (t in 1:(T-1)) {
    log_pi[t+1] <- log_pi[t] + 0.5 * d * (log(T-t) - log(T-t+1))
  }
  return(log_pi - max(log_pi))
}

#' Log Weighted Var-SCP Prior
#'
#' @description
#' Log weighted prior for the Var-SCP model as described in Appendix C.2 of
#' Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that in the
#' absence of a change-point, the posterior probabilities are approximately
#' uniform.
#'
#'
#' @param T An integer. Number of observations in \eqn{y_{1:T}}.
#'
#' @return A numeric vector. Log prior probabilities.
#' @export
#'
log_var_prior <- function(T) {
  log_pi <- rep(0, T)
  for (t in 1:(T-1)) {
    log_pi[t+1] <- sum(log_pi[t] + 0.5,
                       lgamma(0.5 * (T-t+1)) - lgamma(0.5 * (T-t)),
                       0.5 * (T-t) * digamma(0.5 * (T-t)),
                       -0.5 * (T-t+1) * digamma(0.5 * (T-t+1)))
  }
  return(log_pi - max(log_pi))
}

#' Log Weighted MeanVar-SCP Prior
#'
#' @description
#' Log weighted prior for the MeanVar-SCP model as described in Appendix C.2 of
#' Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that in the
#' absence of a change-point, the posterior probabilities are approximately
#' uniform.
#'
#'
#' @param T An integer. Number of observations in \eqn{y_{1:T}}.
#'
#' @return A numeric vector. Log prior probabilities.
#' @export
#'
log_meanvar_prior <- function(T) {
  log_pi <- rep(0, T-1)
  for (t in 1:(T-2)) {
    log_pi[t+1] <- sum(log_pi[t] + 0.5,
                       0.5 * (log(T-t) - log(T-t+1)),
                       lgamma(0.5 * (T-t+1)) - lgamma(0.5 * (T-t)),
                       0.5 * (T-t) * digamma(0.5 * (T-t-1)),
                       -0.5 * (T-t+1) * digamma(0.5 * (T-t)))
  }
  return(c(log_pi - max(log_pi), -Inf))
}
