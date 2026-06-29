#' Multivariate Multiple Independent Change-Point (MICH) Model
#'
#' @description
#' Fits the multivariate version of the MICH model for mean change-points.
#' Number of change-points can either be fixed or `mich_matrix()` will search
#' for the number of changes that maximizes the ELBO when `L_auto == TRUE`.
#'
#' @param y A numeric matrix. T x d matrix of observations.
#' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
#'   intercept is estimated, otherwise `mu_0 = rep(0,d)`.
#' @param fit_scale A logical. If `fit_scale == TRUE`, then the precision matrix
#'   is estimated using the inverse of `var(diff(y))`, otherwise it is assumed
#'   that the precision matrix is equal to `diag(d)`.
#' @param standardize A logical. If `standardize == TRUE`, then `y` is centered
#'    and rescaled before fitting.
#' @param L An integer. Number of mean-scp components included in model. If
#'   `L_auto == TRUE` then `L` lower bounds the number of change-points in the
#'   model.
#' @param L_auto A logical. If `L_auto == TRUE`, then `mich_matrix()` returns
#'   the number of changes between `L` and `L_max` that maximizes the ELBO
#'   (see Appendix C.4 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param L_max L An integer. If `L_auto == TRUE` then `L_max` upper bounds the
#'   number of change-points included in the model.
#' @param pi_l_weighted A logical. If `pi_l_weighted == TRUE`, then the weighted
#'  priors specified in Appendix C.2 of Berlind, Cappello, Madrid Padilla (2025)
#'  are used.
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#' @param verbose A logical. If `verbose == TRUE` and `L_auto == FALSE`, then
#'   the value of the ELBO is printed every 5000th iteration. If
#'   `verbose == TRUE` and `L_auto == TRUE`, the the value of the ELBO is
#'   printed for each L as `mich_matrix()` searches over \[`L`, `L_max`\].
#' @param max_iter An integer. Maximum number of iterations. If ELBO does not
#'   converge before `max_iter` is reached, then `converged == FALSE` in the
#'   returned fit object.
#' @param reverse A logical. If `reverse == TRUE` then MICH is fit to
#'    `y[T:1,]` and the model parameters are reversed in post-processing.
#' @param merge_level A scalar. A value between (0,1) for the significance level
#'   to construct credible sets at when merging. A model component is only
#'   considered to be a candidate for merging if the marginal probability it
#'   represents a change is greater than `merge_level`.
#' @param merge_prob A scalar. A value between (0,1) indicating the merge
#'   criterion. If the posterior probability that two components identify the
#'   same change is greater than `merge_level`, then those components are merged
#'   (see Appendix C.3 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param restart A logical. If `restart == TRUE` and `L_auto == TRUE` then
#'   after `n_search` increments of `L`, if the ELBO has not increased,
#'   `mich_matrix()` will restart by setting the `L` components to the null
#'   model initialization (except for the components with maximum posterior
#'   probabilities > 0.9) then refit and begin the search again.
#' @param n_search An integer. Grid search parameter. Number of times to
#'   increment `L` before terminating automatic procedure when `L_auto == TRUE`.
#' @param increment An integer. Number of components to increment `L` by when
#'   `L_auto == TRUE`.
#' @param omega_l A scalar. Prior precision parameter for mean-scp components of
#'   model.
#' @param log_pi_l A numeric matrix. T x L matrix of prior log change-point
#'   location probabilities for each of the L mean change-points.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution, including:
#'   * `y`: A numeric matrix. Original data.
#'   * `Sigma`: A numeric matrix. Estimate of the precision.
#'   * `L`: An integer. Number of components included in model.
#'   * `pi_bar`: A numeric matrix. A T x L matrix of posterior change-point
#'     location probabilites.
#'   * `residual`: A numeric matrix. Residual `y - mu`.
#'   * `mu`: A numeric matrix. Posterior estimate of mean signal.
#'   * `mu_0`: A numeric vector. Estimate of the intercept.
#'   * `post_params`: A list. List of posterior parameters for each mean-scp
#'     component.
#'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
#'   * `converged`: A logical. Indicates whether relative increase in the ELBO
#'     is less than `tol`.
#'
mich_matrix <- function(
  y, fit_intercept, fit_scale, standardize,
  L, L_auto, L_max, pi_l_weighted,
  tol, verbose, max_iter, reverse,
  merge_level, merge_prob,
  restart, n_search, increment,
  omega_l, log_pi_l
) {
  #### set up ####
  # store original y
  y_raw <- y

  # calculate dimensions of y
  T <- nrow(y)
  d <- ncol(y)

  # log params
  log_omega_l <- log(omega_l)
  d_log_omega_l <- d * log_omega_l

  # min prob to keep component when restarting
  keep_level <- 0.9

  # standardize data for numerical stability
  if (standardize) {
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
  lambda <- 1 / Sigma_eigen$values
  log_lambda <- -log(Sigma_eigen$values)

  # eigen value parameters
  eigen_vals <- omega_l + outer((T-1:(T+1)+1), lambda)
  log_prob_weights <- log_pi_l - 0.5 * rowSums(log(eigen_vals))
  Lambda_bar_log_det <- rowSums(log(eigen_vals))
  inv_Lambda_bar_trace <- omega_l / rowSums(eigen_vals)
  mean_weights <- matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE) / eigen_vals
  sandwich_weights <- mean_weights * matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE)

  # initializing posterior parameters
  post_params <- param_init(L, T, d)

  # initialize mu_0
  if (fit_intercept) {
    mu_0 <- colMeans(y[1:ceiling(log(T)),,drop=FALSE])
  } else mu_0 <- rep(0.0, d)

  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  while (!merged) {
    fit <- mich_matrix_cpp(
      y, mu_0,
      lambda, log_lambda, Q,
      Lambda_bar_log_det, inv_Lambda_bar_trace,
      mean_weights, sandwich_weights, log_prob_weights,
      fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
      omega_l, log_omega_l, d_log_omega_l, log_pi_l,
      post_params
    )

    fit <- multi_component_merge(
      fit, lambda, mean_weights, merge_level, merge_prob
    )

    L <- fit$L
    post_params <- fit$post_params

    if (L_auto) {
      log_pi_l <- log_pi_l[,1:max(1,L), drop=FALSE]
      log_prob_weights <- log_prob_weights[,1:max(1,L), drop=FALSE]
    }

    merged <- fit$merged
    if (verbose & !merged) print(paste0("Merging to L = ", L))
  }

  #### auto procedure ####
  fit <- multi_auto_search(
    y, fit, omega_l, log_pi_l,
    L_max, L_auto, fit_intercept, fit_scale,  max_iter, tol, verbose,
    merge_level, merge_prob, restart, n_search, increment
  )

  #### reverse model if reverse == TRUE ####
  if (reverse) {
    lambda <- fit$lambda
    log_lambda <- log(fit$lambda)
    Q <- fit$Q

    # eigen value parameters
    eigen_vals <- omega_l + outer((T-1:(T+1)+1), lambda)
    log_prob_weights <- log_pi_l - 0.5 * rowSums(log(eigen_vals))
    Lambda_bar_log_det <- rowSums(log(eigen_vals))
    inv_Lambda_bar_trace <- omega_l / rowSums(eigen_vals)
    mean_weights <- matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE) / eigen_vals
    sandwich_weights <- mean_weights * matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE)

    if(verbose) print("reversing model")
    y <- y[T:1,]
    y_raw <- y_raw[T:1,]
    L <- fit$L

    # reversed residuals, variance, and intercepts
    post_params <- fit$post_params
    QTr_bar <- fit$QTr[T:1,]
    mu_0 <- rep(0.0, d)
    for (l in seq_len(L)) {
      mu_0 <- mu_0 + Q %*% post_params[[l]][["QTmu_bar"]][T,]
    }

    # don't reverse weighted priors
    if (!pi_l_weighted) {
      log_pi_l[T:1,] <- log_pi_l[T:1,1:max(1,L), drop = FALSE]
    } else {
      log_pi_l <- sapply(1:max(1,L), function(i) log_pi_l[,1])
    }

    # reversing mean components
    if (L > 0) pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

    for (l in seq_len(L)) {
      tau_l <- which.max(pi_bar_l[,l])
      QTmu_bar_l <- multi_mu_bar_fn(post_params[[l]][["QTb_bar"]], post_params[[l]][["pi_bar"]])
      QTr_bar_l <- QTr_bar + QTmu_bar_l - matrix(colMeans(QTmu_bar_l[tau_l:T,,drop = FALSE]), nrow = T, ncol = d, byrow = TRUE)
      QTr_bar_l <- QTr_bar_l[T:1,]

      post_params[[l]] <- multi_mean_scp(
        QTr_bar_l, lambda, mean_weights, sandwich_weights, log_prob_weights[,l]
      )

      pi_bar_l[,l] <- post_params[[l]][["pi_bar"]]
    }

    # refit model and merge components
    merged <- FALSE # flag indicating components have been merged
    while (!merged) {
      fit <- mich_matrix_cpp(
        y, mu_0,
        lambda, log_lambda, Q,
        Lambda_bar_log_det, inv_Lambda_bar_trace,
        mean_weights, sandwich_weights, log_prob_weights,
        fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
        omega_l, log_omega_l, d_log_omega_l, log_pi_l,
        post_params
      )

      fit <- multi_component_merge(
        fit, lambda, mean_weights, merge_level, merge_prob
      )

      L <- fit$L
      post_params <- fit$post_params

      if (L_auto) {
        log_pi_l <- log_pi_l[,1:max(1,L), drop=FALSE]
        log_prob_weights <- log_prob_weights[,1:max(1,L), drop=FALSE]
      }

      merged <- fit$merged
      if (verbose & !merged) print(paste0("Merging to L = ", L))
    }

    fit <- multi_auto_search(
      y, fit, omega_l, log_pi_l,
      L_max, L_auto, fit_intercept, fit_scale,  max_iter, tol, verbose,
      merge_level, merge_prob, restart, n_search=1, increment
    )
  }

  #### return model ####
  class(fit) <- "mich.fit"
  fit$y <- y_raw

  # calculate correlated parameters
  for (l in seq_len(fit$L)) {
    fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["QTb_bar"]] %*% t(fit$Q)
    fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["QTmu_bar"]] %*% t(fit$Q)
  }
  fit$Sigma <- fit$Q %*% diag(1 / fit$lambda) %*% t(fit$Q)

  # rescale data to original units
  if (standardize) {
    # rescale and center posterior parameters
    fit$mu_0 <- c(fit$mu_0 %*% scale) + center
    fit$mu_bar <- matrix(fit$mu_0, nrow = T, ncol = d, byrow = TRUE)
    for (l in seq_len(fit$L)) {
      fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["b_bar"]] %*% scale
      fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["mu_bar"]] %*% scale
    }
    fit$Sigma <- scale %*% fit$Sigma %*% scale
  }

  if (fit$L > 0) fit$pi_bar <- sapply(1:fit$L, function(l) fit$post_params[[l]][["pi_bar"]])

  # calculate mean and residual
  fit$mu_bar <- matrix(fit$mu_0, nrow = T, ncol = d, byrow = TRUE)
  for (l in seq_len(fit$L)) {
    fit$mu_bar <- fit$mu_bar + fit$post_params[[l]][["mu_bar"]]
  }

  return(fit)
}
