#' Multiple Independent Change-Point (MICH) Model
#'
#' @description
#' Implementation of the MICH model as described in Berlind, Cappello, and
#' Madrid Padilla. MICH is a Bayesian change-point detection method that
#' can quantify uncertainty around estimated change-points in the form of
#' credible sets.
#'
#' @param y A numeric vector or matrix. A length \eqn{T} vector of observations
#'   exhibiting change-points in the mean and/or variance of the series, or a
#'   \eqn{T \times d} matrix of observations exhibiting just mean changes.
#' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
#'   intercept is estimated, otherwise it is assumed  that
#'   \eqn{\mu_0 = 0}.
#' @param fit_scale A logical. If `fit_scale == TRUE`, then the initial
#'   precision is estimated, otherwise it is assumed  that \eqn{\lambda_0 = 1}
#'   (or \eqn{\Lambda_0 = \mathbf{I}_d} if `y` is a matrix).
#' @param standardize A logical. If `standardize == TRUE`, then `y` is centered
#'    and rescaled before fitting.
#' @param L,J,K Integers. Respective number of mean-scp, var-scp, and
#'   meanvar-scp components included in the model. If `L_auto == TRUE`,
#'   `K_auto == TRUE`, or `J_auto == TRUE` then `L`, `K`, and `J` lower bound
#'   the number of each kind of change-point in the model.
#' @param L_auto,K_auto,J_auto Logicals. If `L_auto == TRUE`,
#'   `K_auto == TRUE`, and/or `J_auto == TRUE`, then `mich_vector()` searches
#'   forr and returns the \eqn{L} between `L` and `L_max`, the \eqn{K} between
#'   `K` and `K_max`, and/or the \eqn{J} between `J` and `J_max` that maximize
#'   the ELBO (see Appendix C.4 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param L_max,K_max,J_max Integers. If `L_auto == TRUE`,
#'   `K_auto == TRUE`, or `J_auto == TRUE` then `L_max`, `K_max`, and `J_max`
#'   upper bound the number of each kind of change-point in the model.
#' @param alpha A scalar. A value in (0,1) controlling the significance level.
#' @param null_prob A scalar. Prior probability that each component does not
#'   contain a change.
#' @param prior A string. One of `c('var.weighted', 'weighted', 'independent')`
#'   specifying the model for the prior precision.
#' @param log_omega A scalar. Prior precision parameter.
#' @param pi_l Prior for the mean-scp model components. Must either be a
#'   length \eqn{T} vector with entries that sum to one, a \eqn{T \times L}
#'   matrix with columns that sum to one, or `NULL`, in which case
#'   `pi_j = c(rep((1-null_prob) / T, T), null_prob)`
#' @param pi_k Prior for the var-scp model components. Must either be a
#'   length \eqn{T} vector with entries that sum to one, a \eqn{T \times K}
#'   matrix with columns that sum to one, or `NULL`, in which case
#'   `pi_j = c(rep((1-null_prob) / T, T), null_prob)`
#' @param pi_j Prior for the meanvar-scp model components. Must either be a
#'   length \eqn{T} vector with entries that sum to one, a \eqn{T \times J}
#'   matrix with columns that sum to one, or `NULL`, in which case
#'   `pi_j = c(rep((1-null_prob) / T, T), null_prob)`
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#'   Once \eqn{(\text{ELBO}_{t+1}-\text{ELBO}_t)/ \text{ELBO}_t} falls below
#'   `tol` the variational algorithm terminates.
#' @param merge_prob A scalar. A value between (0,1) indicating the merge
#'   criterion. If the posterior probability that two components identify the
#'   same change is greater than `merge_prob`, then those components are merged.
#' @param max_iter An integer. Maximum number of variational iterations. If ELBO
#'   does not converge before `max_iter` is reached, then `converged == FALSE`
#'   in the returned fit object.
#' @param verbose A logical. If `verbose == TRUE` and `L_auto == FALSE`,
#'   `K_auto == FALSE`, and `J_auto == FALSE` then the value of the ELBO is
#'   printed every 5000th iteration. If `verbose == TRUE` and any of `L_auto`,
#'   `K_auto`, or `J_auto` are `TRUE`, then the value of the ELBO is printed for
#'   each combinatin of \eqn{(L,K,J)} as `mich_vector()` searches for the
#'   parameterization that maximized the EBLO.
#' @param reverse A logical. If `reverse == TRUE` then MICH is fit to
#'   \eqn{\mathbf{y}_{T:1}} and the model parameters are reversed in
#'   post-processing.
#' @param restart logical. If `restart == TRUE` and `L_auto`, `K_auto`, or
#'   `J_auto` are `TRUE`, then after `n_search` increments of `L`, `K`, and/or
#'   `J`, if the ELBO has not increased, `mich()` will restart by setting
#'   the `L`, `K`, and `J` components to the null model initialization then
#'   refit and begin the search again.
#' @param increment An integer. Number of components to increment `L`, `K` and
#'   `J` by when `L_auto`, `K_auto`, or `J_auto` are `TRUE`.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution. If `y` is a vector, this list includes:
#'   * `y`: A numeric vector. Original data.
#'   * `L`,`K`,`J`: Integers. Number of mean-scp, var-scp, and meanvar-scp
#'     components included in the model.
#'   * `residual`: A numeric vector. Residual \eqn{\tilde{\mathbf{r}}_{1:T}}
#'     (see (B.4) of Berlind, Cappello, Madrid Padilla (2025)).
#'   * `mu`: A numeric vector. Posterior estimate of
#'     \eqn{\Sigma_{\ell=1}^L E[\mu_{\ell,1:T}|\mathbf{y}_{1:T}] + \Sigma_{j=1}^J E[\mu_{j,1:T}|\mathbf{y}_{1:T}]}.
#'   * `lambda`: A numeric vector. Posterior estimate of
#'     \eqn{\Pi_{k=1}^K E[\lambda_{k,1:T}|\mathbf{y}_{1:T}] \times \Pi_{j=1}^J E[\lambda_{j,1:T}|\mathbf{y}_{1:T}]}.
#'   * `delta`: A numeric vector. Posterior estimate of equation (B.4) of
#'     Berlind, Cappello, Madrid Padilla (2025).
#'   * `mu_0`: A scalar. Estimate of the intercept.
#'   * `lambda_0`: A scalar. Estimate of the initial precision.
#'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
#'   * `converged`: A logical. Indicates whether relative increase in the ELBO
#'     is less than `tol`.
#'   * `meanvar_model`: A list. List of meanvar-scp posterior parameters:
#'     * `pi_bar`: A numeric matrix. A \eqn{T \times J} matrix of posterior
#'       change-point location probabilities.
#'     * `b_bar`: A numeric matrix. A \eqn{T \times J} matrix of posterior
#'       mean parameters.
#'     * `omega_bar`: A numeric matrix. A \eqn{T \times J} matrix of posterior
#'       precision parameters.
#'     * `v_bar`: A numeric matrix. A \eqn{T \times J} matrix of posterior
#'       rate parameters.
#'     * `u_bar`: A numeric vector. A length \eqn{T} vector of posterior shape
#'       parameters.
#'     * `mu_lambda_bar`: A numeric matrix. A \eqn{T \times J} matrix of scaled
#'       posterior mean signals \eqn{E[\lambda_{jt}\mu_{jt}|\mathbf{y}_{1:T}]}.
#'     * `mu2_lambda_bar`: A numeric matrix. A \eqn{T \times J} matrix of scaled
#'       posterior squared mean signals
#'       \eqn{E[\lambda_{jt}\mu^2_{jt}|\mathbf{y}_{1:T}]}.
#'     * `lambda_bar`: A numeric matrix. A \eqn{T \times J} matrix of
#'       posterior precision signals \eqn{E[\lambda_{jt}|\mathbf{y}_{1:T}]}.
#'   * `mean_model`: A list. List of mean-scp posterior parameters:
#'     * `pi_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior
#'       change-point location probabilities.
#'     * `b_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior
#'       mean parameters.
#'     * `omega_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior
#'       precision parameters.
#'     * `mu_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior mean
#'       signals \eqn{E[\mu_{\ell t}|\mathbf{y}_{1:T}]}.
#'     * `mu2_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior
#'       squared mean signals \eqn{E[\mu_{\ell t}^2|\mathbf{y}_{1:T}]}.
#'   * `var_model`: A list. List of var-scp posterior parameters:
#'     * `pi_bar`: A numeric matrix. A \eqn{T \times K} matrix of posterior
#'       change-point location probabilities.
#'     * `v_bar`: A numeric matrix. A \eqn{T \times K} matrix of posterior
#'       rate parameters.
#'     * `u_bar`: A numeric vector. A length \eqn{T} vector of posterior shape
#'       parameters.
#'     * `lambda_bar`: A numeric matrix. A \eqn{T \times K} matrix of
#'       posterior precision signals \eqn{E[\lambda_{kt}|\mathbf{y}_{1:T}]}.
#'
#' If `y` is a matrix, this list includes:
#'   * `y`: A numeric matrix. Original data.
#'   * `Sigma`: A numeric matrix. Estimate of \eqn{\text{var}(y_t)}.
#'   * `lambda`: A numeric vector. Common eigenvalues of \eqn{\text{Var}(y_t)}.
#'   * `Q`: A numeric matrix. Common eigenvectors of \eqn{\text{Var}(y_t)}.
#'   * `L`: An integer. Number of components included in model.
#'   * `pi_bar`: A numeric matrix. A \eqn{T \times L} matrix of posterior change-point
#'     location probabilites.
#'   * `QTr`: A numeric matrix. Decorrelated residual \eqn{Q^T(y_t - \mu_t)}.
#'   * `mu_bar`: A numeric matrix. Posterior estimate of mean signal.
#'   * `mu_0`: A numeric vector. Estimate of the intercept.
#'   * `post_params`: A list. List of posterior parameters for each mean-scp
#'     component.
#'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
#'   * `converged`: A logical. Indicates whether relative increase in the ELBO
#'     is less than `tol`.
#'   * `merged`: A logical. Flag indicating whether any components were merged.
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' }
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
#' fit = mich(Y, L_auto = TRUE, null_prob = 0.05, log_omega = -5)
#' plot(fit, level = 0.95, cs = TRUE, signal = TRUE)
#'
mich <- function(
    y, fit_intercept = TRUE, fit_scale = TRUE, standardize = TRUE,
    J = 0, L = 0, K = 0,
    J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
    J_max = Inf, L_max = Inf, K_max = Inf,
    alpha = 0.95, null_prob = 1-alpha,
    prior = c('var.weighted', 'weighted', 'independent'),
    log_omega = NULL, pi_l=NULL, pi_k=NULL, pi_j=NULL,
    tol = 1e-5, merge_prob = NULL,
    max_iter = 1e4, verbose = FALSE, reverse = FALSE,
    restart = TRUE, increment = 1
) {
  #### data validation ####
  # checking that y is a numeric vector/matrix
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.numeric(y)) stop("y must contain only numeric data.")
  if (!is.vector(y) & !is.matrix(y)) stop("y must be a numeric vector or matrix.")

  # calculate dimensions of y and handle NAs
  if (all(is.na(y)) | length(y) == 0) stop("y does not contain any data.")
  if (is.matrix(y)) {
    if (any(is.na(y))) {
      warning("NAs detected in y. Removing missing rows and columns and replacing missing values with preceding observation.")
      y <- y[!apply(y, 1, function(x) all(is.na(x))), !apply(y, 2, function(x) all(is.na(x))), drop == FALSE]
      for (t in 1:nrow(y)) {
        for (i in 1:ncol(y)) {
          if (is.na(y[t, i])){
            if (t == 1) y[t, i] <- y[which.min(is.na(y[,i])),i]
            else y[t, i] <- y[t-1, i]
          }
        }
      }
    }
    d <- ncol(y)
    T <- nrow(y)
    if (d == 1) y <- as.vector(y)
  } else {
    if (any(is.na(y))) {
      warning("NAs detected in y. Removing missing values.")
      y <- y[!is.na(y)]
    }
    d <- 1
    T <- length(y)
  }

  if (T < 2) stop("y must have at least 2 observations.")

  #### merge and detect defaults ####
  alpha <- scalar_check(alpha)
  if (alpha < 0 | alpha > 1) stop("alpha must be in [0,1].")
  if (is.null(merge_prob)) merge_prob <- (log(T) / T)^2
  merge_prob <- scalar_check(merge_prob)
  if (merge_prob < 0 | merge_prob > 1) stop("merge_prob must be in [0,1].")
  n_search <- max(4, ceiling(log(T) / ((1 + restart) * increment)))

  #### checking other model parameters ####
  prior <- match.arg(prior)
  max_iter <- integer_check(max_iter)
  tol <- scalar_check(tol)
  fit_intercept <- logical_check(fit_intercept)
  fit_scale <- logical_check(fit_scale)
  verbose <- logical_check(verbose)

  J_auto <- logical_check(J_auto)
  L_auto <- logical_check(L_auto)
  K_auto <- logical_check(K_auto)

  J <- integer_check(J)
  L <- integer_check(L)
  K <- integer_check(K)

  J_max <- integer_check(J_max)
  L_max <- integer_check(L_max)
  K_max <- integer_check(K_max)

  if (J_max < J) "J_max must be >= J."
  if (L_max < L) "L_max must be >= L."
  if (K_max < K) "K_max must be >= K."

  if (!J_auto) J_max <- J
  if (!L_auto) L_max <- L
  if (!K_auto) K_max <- K

  if (J > 0 & J_auto & verbose) warning(paste0("J_auto = TRUE and J > 0. Beginning search from J = ", J))
  if (L > 0 & L_auto & verbose) warning(paste0("L_auto = TRUE and L > 0. Beginning search from L = ", L))
  if (K > 0 & K_auto & verbose) warning(paste0("K_auto = TRUE and K > 0. Beginning search from K = ", K))

  ##### checking that priors are proper ####
  null_prob <- scalar_check(null_prob)
  if (null_prob < 0 | null_prob > 1) stop("null_prob must be in [0,1].")
  if (!is.null(log_omega)) {
    log_omega <- scalar_check(log_omega)
    omega <- exp(log_omega)
    #if(verbose) warning(paste0("omega = ", omega, ". This may not result in Type-I error control at alpha = ", alpha))
    log_omega_j <- log_omega; u_j <- omega; v_j <- omega;
    log_omega_l <- log_omega;
    u_k <- omega; v_k <- omega;
  }

  # J components
  if (!is.null(pi_j)) {
    pi_j <- prob_check(pi_j, J, T + 1)
    if(verbose) warning(paste0("Ignoring null_prob and using provided pi_j."))
  } else {
    pi_j <- matrix((1-null_prob) / T, nrow = T, ncol = max(1,J))
    pi_j <- rbind(pi_j, null_prob)
  }
  log_pi_j <- log(pi_j)

  # L components
  if (!is.null(pi_l)) {
    pi_l <- prob_check(pi_l, L, T + 1)
    if(verbose) warning(paste0("Ignoring null_prob and using provided pi_l."))
  } else {
    pi_l <- matrix((1-null_prob) / T, nrow = T, ncol = max(1,L))
    pi_l <- rbind(pi_l, null_prob)
  }
  log_pi_l <- log(pi_l)

  # K components
  if (!is.null(pi_k)) {
    pi_k <- prob_check(pi_k, K, T)
    if(verbose) warning(paste0("Ignoring null_prob and using provided pi_k."))
  } else {
    pi_k <- matrix((1-null_prob) / T, nrow = T, ncol = max(1,K))
    pi_k <- rbind(pi_k, null_prob)
  }
  log_pi_k <- log(pi_k)

  #### fit model ####
  if (is.matrix(y)) {
    # mulitvariate model ####
    if (J > 0 | K > 0) Warning("MICH currently only suports mean changepoint detection for multivariate y. Setting J = K = 0.")
    if (T < d & fit_scale) {
      warning("y has more columns than rows. MICH does not currently support high-dimensional variance estimation. Assuming Var(y) = I.")
      fit_scale <- FALSE
    }

    # standardize data
    if (standardize){
      fit_intercept <- TRUE
      y_raw <- y

      center <- colMeans(y)
      scale_eigen <- eigen(var(y))
      Q_scale <- scale_eigen$vectors
      lambda_scale <- scale_eigen$values

      if (any(lambda_scale <= 0)) {
        warning("Var(y) is singular. Consider removing collinear columns.")
        lambda_scale <- lambda_scale + 1e-5
      }

      scale <- Q_scale %*% diag(sqrt(lambda_scale)) %*% t(Q_scale)
      inv_scale <- Q_scale %*% diag(1 / sqrt(lambda_scale)) %*% t(Q_scale)

      y <- (y - matrix(center, nrow = T, ncol = d, byrow = TRUE)) %*% inv_scale
    }

    if (fit_scale) {
      # estimate precision matrix
      y_diff <- diff(y) # difference out mean changes

      # remove outliers due to big mean changes
      y_diff_norm <- sqrt(rowSums(y_diff^2))
      y_diff <- y_diff[y_diff_norm <= stats::quantile(y_diff_norm, p = 0.75) +  1.5 * stats::IQR(y_diff_norm), ]

      # estimate variance
      Sigma <- (t(y_diff) %*% y_diff) / (2 * nrow(y_diff))
      Sigma_eigen <- eigen(Sigma)
      if (any(Sigma_eigen$values <= 0)) {
        warning("Var(y) is singular. Consider removing collinear columns.")
        Sigma_eigen$values <- Sigma_eigen$values + 1e-5
      }
      Q <- Sigma_eigen$vectors
      lambda <- 1.0 / Sigma_eigen$values
    } else {
      Q <- diag(d)
      lambda <- rep(1.0, d)
    }

    # determine omega_l
    if (is.null(log_omega)) {
      log_omega_l <- log_omega_mean_fdr(alpha, N = 1000, T, prior, log_pi_l[,1], lambda)
    }

    # fit model
    if (!reverse) {
      fit <- mich_matrix(
        y, lambda, Q,
        fit_intercept, fit_scale,
        L, L_auto, L_max, prior,
        tol, verbose, max_iter,
        merge_prob, restart, n_search, increment, reverse,
        log_omega_l, log_pi_l
      )
    } else {
      log_pi_l[1:T,] <- log_pi_l[T:1,] # reverse prior
      fit <- mich_matrix(
        y[T:1,], lambda, Q,
        fit_intercept=TRUE, fit_scale,
        L, L_auto, L_max, prior,
        tol, verbose, max_iter,
        merge_prob, restart, n_search, increment, reverse,
        log_omega_l, log_pi_l
      )
    }

    # rescale data to original units
    if (standardize) {
      # rescale and center posterior parameters
      fit$y <- y_raw
      fit$mu_0 <- c(fit$mu_0 %*% scale) + center
      fit$mu_bar <- fit$mu_bar %*% scale +  matrix(center, nrow = T, ncol = d, byrow = TRUE)
      for (l in seq_len(fit$L)) {
        fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["b_bar"]] %*% scale
        fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["mu_bar"]] %*% scale
      }
      fit$Sigma <- scale %*% fit$Sigma %*% scale
    }
  } else {
    if (reverse) {
      if (!pi_l_weighted) log_pi_l[1:T,] <- log_pi_l[T:1,,drop = FALSE]
      if (!pi_k_weighted) log_pi_k[1:T,] <- log_pi_k[T:1,,drop = FALSE]
      if (!pi_j_weighted) log_pi_j[1:T,] <- log_pi_j[T:1,,drop = FALSE]

      mich_vector(
        y[T:1], fit_intercept = TRUE, fit_scale = TRUE, standardize,
        J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
        pi_j_weighted, pi_l_weighted, pi_k_weighted,
        tol, verbose, max_iter, reverse,
        merge_level, merge_prob,
        restart, n_search, increment,
        omega_j, u_j, v_j, log_pi_j,
        omega_l, log_pi_l,
        u_k, v_k, log_pi_k
      )
    } else {
      mich_vector(
        y, fit_intercept, fit_scale, standardize,
        J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
        pi_j_weighted, pi_l_weighted, pi_k_weighted,
        tol, verbose, max_iter, reverse,
        merge_level, merge_prob,
        restart, n_search, increment,
        omega_j, u_j, v_j, log_pi_j,
        omega_l, log_pi_l,
        u_k, v_k, log_pi_k
      )
    }
  }
  return(fit)
}

