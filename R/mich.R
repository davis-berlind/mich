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
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#'   Once \eqn{(\text{ELBO}_{t+1}-\text{ELBO}_t)/ \text{ELBO}_t} falls below
#'   `tol` the variational algorithm terminates.
#' @param merge_prob A scalar. A value between (0,1) indicating the merge
#'   criterion. If the posterior probability that two components identify the
#'   same change is greater than `merge_level`, then those components are merged
#'   (see Appendix C.3 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param merge_level A scalar. A value between (0,1) for the significance level
#'   to construct credible sets at when merging. A model component is only
#'   considered to be a candidate for merging if its `merge_level`-level
#'   credible set contains fewer than `detect` indices.
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
#'   `J`, if the ELBO has not increased, `mich_vector()` will restart by setting
#'   the `L`, `K`, and `J` components to the null model initialization (except
#'   for the components with maximum posterior probabilities > 0.9) then refit
#'   and begin the search again.
#' @param increment An integer. Number of components to increment `L`, `K` and
#'   `J` by when `L_auto`, `K_auto`, or `J_auto` are `TRUE`.
#' @param omega_j,u_j,v_j Scalars. Prior precision, shape, and rate parameters
#'   for meanvar-scp components of the model.
#' @param pi_j A character or numeric vector or matrix. Prior for the
#'   meanvar-scp model components. Must either be a length \eqn{T} vector with
#'   entries that sum to one, a \eqn{T \times J} matrix with columns that sum to
#'   one, or a character equal to `"weighted"` or `"uniform"`. If `pi_j` is a
#'   vector and `J > 1`, then `pi_j` is used as the prior for all `J`
#'   meanvar-scp components in the model. If `pi_j == "uniform"` then the
#'   uniform prior \eqn{\pi_{jt} = 1/T} is used for all components. If
#'   `pi_j == "weighted"` then `log_meanvar_prior()` is used to calculate the
#'   weighted prior as described in Appendix C.2 of Berlind, Cappello, Madrid
#'   Padilla (2025). If `J_auto == TRUE` then it must be the case that
#'   `pi_j %in% c("uniform", "weighted")`.
#' @param omega_l A scalar. Prior precision parameter for mean-scp components of
#'   the model. If `y` is a matrix then the prior precision is
#'   \eqn{\omega_\ell\mathbf{I}_d}.
#' @param pi_l A character or numeric vector or matrix. Prior for the mean-scp
#'   model components. Must either be a length \eqn{T} vector with entries that
#'   sum to one, a \eqn{T \times L} matrix with columns that sum to one, or a
#'   character equal to `"weighted"` or `"uniform"`. If `pi_l` is a vector and
#'   `L > 1`, then `pi_l` is used as the prior for all `L` mean-scp components.
#'   If `pi_l == "uniform"` then the uniform prior \eqn{\pi_{\ell t} = 1/T} is
#'   used for all components. If `pi_l == "weighted"` then `log_mean_prior()`
#'   is used to calculate the weighted prior as described in Appendix C.2 of
#'   Berlind, Cappello, Madrid Padilla (2025). If `L_auto == TRUE` then it must
#'   be the case that `pi_l %in% c("uniform", "weighted")`.
#' @param u_k,v_k Scalar. Prior shape and rate parameters for var-scp components
#'   of the model.
#' @param pi_k A character or numeric vector or matrix. Prior for the var-scp
#'   model components. Must either be a length \eqn{T} vector with entries that
#'   sum to one, a \eqn{T \times K} matrix with columns that sum to one, or a
#'   character equal to `"weighted"` or `"uniform"`. If `pi_k` is a vector and
#'   `K > 1`, then `pi_k` is used as the prior for all `K` var-scp components.
#'   If `pi_k == "uniform"` then the uniform prior \eqn{\pi_{k t} = 1/T} is used
#'   for all components. If `pi_k == "weighted"` then `log_var_prior()` is used
#'   to calculate the weighted prior as described in Appendix C.2 of Berlind,
#'   Cappello, Madrid Padilla (2025). If `K_auto == TRUE` then it must be the
#'   case that `pi_k %in% c("uniform", "weighted")`.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution. If `y` is a vector, this list includes:
#'   * `y`: A numeric vector. Original data.
#'   * `L`,`K`,`J`: Integers. Number of mean-scp, var-scp, and meanvar-scp
#'     components included in the model.
#'   * `residual`: A numeric vector. Residual \eqn{\tilde{\mathbf{r}}_{1:T}}
#'     (see (B.4) of Berlind, Cappello, Madrid Padilla (2025)).
#'   * `mu`: A numeric vector. Posterior estimate of
#'     \eqn{\Sigma_{\ell=1}^L E[\boldsymbol{\mu}_{\ell,1:T}|\mathbf{y}_{1:T}] + \Sigma_{j=1}^J E[\boldsymbol{\mu}_{j,1:T}|\mathbf{y}_{1:T}]}.
#'   * `lambda`: A numeric vector. Posterior estimate of
#'     \eqn{\Pi_{k=1}^K E[\boldsymbol{\lambda}_{k,1:T}|\mathbf{y}_{1:T}] \times \Pi_{j=1}^J E[\boldsymbol{\lambda}_{j,1:T}|\mathbf{y}_{1:T}]}.
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
#' If `y` is a vector, this list includes:
#'   * `y`: A numeric matrix. Original data.
#'   * `Sigma`: A numeric matrix. Estimate of \eqn{\Lambda^{-1}} if
#'     `fit_scale == TRUE`.
#'   * `L`: An integer. Number of components included in model.
#'   * `pi_bar`: A numeric matrix. A  \eqn{T \times L} matrix of posterior
#'     change-point location probabilites.
#'   * `residual`: A numeric matrix. Residual \eqn{\mathbf{r}_{1:T}} after
#'     subtracting out each \eqn{E[\boldsymbol{\mu}_{\ell t}]} from
#'     \eqn{\mathbf{y}_{1:T}}.
#'   * `mu`: A numeric matrix. Posterior estimate of
#'     \eqn{\Sigma_{\ell=1}^L E[\boldsymbol{\mu}_{\ell,1:T}|\mathbf{y}_{1:T}]}.
#'   * `mu_0`: A numeric vector. Estimate of the intercept.
#'   * `post_params`: A list. List of posterior parameters for each mean-scp
#'     component.
#'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
#'   * `converged`: A logical. Indicates whether relative increase in the ELBO
#'     is less than `tol`.
#'
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
mich <- function(y, fit_intercept = TRUE, fit_scale = TRUE, standardize = TRUE,
                 J = 0, L = 0, K = 0,
                 J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                 J_max = Inf, L_max = Inf, K_max = Inf,
                 tol = 1e-5, merge_prob = NULL, merge_level = 0.95,
                 max_iter = 1e4, verbose = FALSE, reverse = FALSE,
                 restart = TRUE, increment = 1,
                 omega_j = 1e-3, u_j = 1e-3, v_j = 1e-3, pi_j = "weighted",
                 omega_l = 1e-3, pi_l = "weighted",
                 u_k = 1e-3, v_k = 1e-3, pi_k = "weighted") {

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
  delta <- 0.5
  merge_level <- scalar_check(merge_level)
  if (merge_level < 0 | merge_level > 1) stop("merge_level must be in [0,1].")
  detect <- ceiling(log(T)^(1 + delta))
  if (is.null(merge_prob)) merge_prob <- detect / T^2
  merge_prob <- scalar_check(merge_prob)
  n_search <- max(4, ceiling(log(T) / ((1 + restart) * increment)))

  #### checking other model parameters ####
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

  if (J > 0 & J_auto) warning(paste0("J_auto = TRUE and J > 0. Beginning search from J = ", J))
  if (L > 0 & L_auto) warning(paste0("L_auto = TRUE and L > 0. Beginning search from L = ", L))
  if (K > 0 & K_auto) warning(paste0("K_auto = TRUE and K > 0. Beginning search from K = ", K))

  ##### checking that priors are proper ####

  # J components
  omega_j <- scalar_check(omega_j)
  u_j <- scalar_check(u_j)
  v_j <- scalar_check(v_j)

  pi_j <- prob_check(pi_j, J, T)
  pi_j_weighted <- (pi_j == "weighted")
  if (J_auto & length(pi_j) > 1) stop("When J_auto = TRUE, pi_j must one of 'uniform' or 'weighted'.")
  if (is.character(pi_j)) {
    if (pi_j == "weighted") log_pi_j <- log_meanvar_prior(T)
    else if (pi_j == "uniform") log_pi_j <- rep(0, T)
    log_pi_j <- sapply(1:max(1, J), function(i) log_pi_j)
  } else log_pi_j <- log(pi_j)

  # L components
  omega_l <- scalar_check(omega_l)

  pi_l <- prob_check(pi_l, L, T)
  pi_l_weighted <- (pi_l == "weighted")
  if (L_auto & length(pi_l) > 1) stop("When L_auto = TRUE, pi_l must one of 'uniform' or 'weighted'.")
  if (is.character(pi_l)) {
    if (pi_l == "weighted") log_pi_l <- log_mean_prior(T, floor(1.5 * d))
    else if (pi_l == "uniform") log_pi_l <- rep(0, T)
    log_pi_l <- sapply(1:max(1, L), function(i) log_pi_l)
  } else log_pi_l <- log(pi_l)

  # K components
  u_k <- scalar_check(u_k)
  v_k <- scalar_check(v_k)

  pi_k <- prob_check(pi_k, K, T)
  pi_k_weighted <- (pi_k == "weighted")
  if (K_auto & length(pi_k) > 1) stop("When K_auto = TRUE, pi_k must one of 'uniform' or 'weighted'.")
  if (is.character(pi_k)) {
    if (pi_k == "weighted") log_pi_k <- log_var_prior(T)
    else if (pi_k == "uniform") log_pi_k <- rep(0, T)
    log_pi_k <- sapply(1:max(1, K), function(i) log_pi_k)
  } else log_pi_k <- log(pi_k)

  # call multivariate or univariate MICH
  if (is.matrix(y)) {
    if (J > 0 | K > 0) stop("MICH currently only suports mean changepoint detection for multivariate y. Try setting J = K = 0.")
    if (T < d & fit_scale) {
      warning("y has more columns than rows. MICH does not currently support high-dimensional variance estimation. Assuming Var(y) = I.")
      fit_scale <- FALSE
    }

    if (reverse) {
      if (!pi_l_weighted) pi_l <- pi_l[T:1,,drop = FALSE]

      mich_matrix(y[T:1,], fit_intercept = TRUE, fit_scale, standardize,
                  L, L_auto, L_max, pi_l_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob,
                  restart, n_search, increment,
                  omega_l, log_pi_l)
    } else {
      mich_matrix(y, fit_intercept, fit_scale, standardize,
                  L, L_auto, L_max, pi_l_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob,
                  restart, n_search, increment,
                  omega_l, log_pi_l)
    }
  } else {
    if (reverse) {
      if (!pi_l_weighted) log_pi_l <- log_pi_l[T:1,,drop = FALSE]
      if (!pi_k_weighted) log_pi_k <- log_pi_k[T:1,,drop = FALSE]
      if (!pi_j_weighted) log_pi_j <- log_pi_j[T:1,,drop = FALSE]

      mich_vector(y[T:1], fit_intercept = TRUE, fit_scale = TRUE, standardize,
                  J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                  pi_j_weighted, pi_l_weighted, pi_k_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob,
                  restart, n_search, increment,
                  omega_j, u_j, v_j, log_pi_j,
                  omega_l, log_pi_l,
                  u_k, v_k, log_pi_k)
    } else {
      mich_vector(y, fit_intercept, fit_scale, standardize,
                  J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                  pi_j_weighted, pi_l_weighted, pi_k_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob,
                  restart, n_search, increment,
                  omega_j, u_j, v_j, log_pi_j,
                  omega_l, log_pi_l,
                  u_k, v_k, log_pi_k)
    }
  }
}

