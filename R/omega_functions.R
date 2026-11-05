log_omega_mean_fdr <- function(alpha = 0.95, N = 1000, T, prior, log_pi, lambda) {
  set.seed(422)
  d <- length(lambda)
  std.dev <- 1 / sqrt(lambda)
  Ys <- lapply(1:N, function(i) t(sapply(1:T, function(i) rnorm(d, sd = std.dev))))

  log_omega <- -1
  tnr <- mean_fdr(log_omega, 0.5, Ys, prior, log_pi, lambda, d, T)
  while(tnr < alpha) {
    log_omega <- log_omega - 1
    tnr <- mean_fdr(log_omega, 0.5, Ys, prior, log_pi, lambda, d, T)
  }
  while(tnr > alpha) {
    log_omega <- log_omega + 1
    tnr <- mean_fdr(log_omega, 0.5, Ys, prior, log_pi, lambda, d, T)
  }
  log_omegas <- seq(log_omega - 1, log_omega, length.out = 10)
  while (TRUE) {
    tnrs <- sapply(log_omegas, function(log_omega) {
      mean_fdr(log_omega, 0.5, Ys, prior, log_pi, lambda, d, T)
    })
    errors <- abs(tnrs - alpha)
    min_dx <- which.min(errors)
    if (errors[min_dx] < 1e-5) return(log_omegas[min_dx])
    log_omegas <- seq(log_omegas[min_dx-1], log_omegas[min_dx+1], length.out=10)
  }
}

mean_fdr <- function(log_omega, alpha, Ys, prior, log_pi, lambda, d, T) {
  # if (is.null(log_pi)) {
  #   logit <- 0.5 * d * (log(omega) - log(1 + omega) +  1 / (1 + omega)) + max(log(T), sqrt(d * log(T)))
  #   log_pi_T1 <- logit - log(-log(1-alpha) + exp(logit))
  #   log_pi <- c(rep(log(1 - exp(log_pi_T1)) - log(T), T), log_pi_T1)
  # }
  if (prior == "var.weighted") {
    omega <- exp(log_omega)
    omegas <- outer(T:1, omega * lambda)
    log_det_omega_omega_bar <- rep(d * (log_omega - log(1+omega)), T)
  } else if (prior == 'weighted') {
    omegas <- omega * matrix(T:1, nrow = T, ncol = d, byrow = FALSE)
    log_det_omega_omega_bar <- rep(d * log_omega - sum(log(lambda + omega)), T)
  } else if (prior == 'independent') {
    omegas <- matrix(omega, nrow = T, ncol = d)
    log_det_omega_omega_bar <- d * log_omega - rowSums(log(outer(T:1, lambda) + omegas))
  }
  eigen_vals <- outer(T:1, lambda) + omegas
  sandwich_weights <- 0.5 * matrix(lambda^2, nrow = T, ncol = d, byrow = TRUE) / eigen_vals
  log_prob <- log_pi + c(log_det_omega_omega_bar, 0.0)

  probs <- unlist(lapply(Ys, function(y) multi_mean_prob(y, lambda, sandwich_weights, log_prob)[T+1]))
  return(mean(probs > (1-alpha)))
}

