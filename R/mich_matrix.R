#' Multivariate Multiple Independent Change-Point (MICH) Model
#'
#' @description
#' Fits the multivariate version of the MICH model for mean change-points.
#' Number of change-points can either be fixed or `mich_matrix()` will search
#' for the number of changes that maximizes the ELBO when `L_auto == TRUE`.
#'
#' @param y A numeric matrix. \eqn{T \times d} matrix of observations.
#' @param lambda A numeric vector. Common eigenvalues of \eqn{\text{Var}(y_t)}.
#' @param Q A numeric matrix. Common eigenvectors of \eqn{\text{Var}(y_t)}.
#' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
#'   intercept is estimated, otherwise `mu_0 = rep(0,d)`.
#' @param fit_scale A logical. If `fit_scale == TRUE`, then the precision matrix
#'   is estimated using the inverse of `var(diff(y))/2`, otherwise it is assumed
#'   that the precision matrix is equal to `diag(d)`.
#' @param L An integer. Number of mean-scp components included in model. If
#'   `L_auto == TRUE` then `L` lower bounds the number of change-points in the
#'   model.
#' @param L_auto A logical. If `L_auto == TRUE`, then `mich_matrix()` returns
#'   the number of changes between `L` and `L_max` that maximizes the ELBO. See
#'   `multi_auto_search()`.
#' @param L_max An integer. If `L_auto == TRUE` then `L_max` upper bounds the
#'   number of change-points included in the model.
#' @param prior A string. One of `c('var.weighted', 'weighted', 'independent')`
#'   specifying the model for the prior precision.
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#' @param verbose A logical. If `verbose == TRUE` and `L_auto == FALSE`, then
#'   the value of the ELBO is printed every 5000th iteration. If
#'   `verbose == TRUE` and `L_auto == TRUE`, the the value of the ELBO is
#'   printed for each `L` as `mich_matrix()` searches over \[`L`, `L_max`\].
#' @param max_iter An integer. Maximum number of iterations. If ELBO does not
#'   converge before `max_iter` is reached, then `converged == FALSE` in the
#'   returned fit object.
#' @param merge_prob A scalar. A value between (0,1) indicating the merge
#'   criterion. If the posterior probability that two components identify the
#'   same change is greater than `merge_prob`, then those components are merged.
#' @param restart A logical. If `restart == TRUE` and `L_auto == TRUE` then
#'   after `n_search` increments of `L`, if the ELBO has not increased,
#'   `mich_matrix()` will restart by setting the `L` components to the null
#'   model initialization then refit and begin the search again.
#' @param n_search An integer. Grid search parameter. Number of times to
#'   increment `L` before terminating automatic procedure when `L_auto == TRUE`.
#' @param increment An integer. Number of components to increment `L` by when
#'   `L_auto == TRUE`.
#' @param reverse A logical. If `reverse == TRUE` then MICH is fit to
#'    `y[T1,]` and the model parameters are reversed in post-processing.
#' @param log_omega_l A scalar. Prior precision parameter.
#' @param log_pi_l A numeric matrix. \eqn{(T+1) \times L} matrix of prior log
#'   change-point location probabilities for each of the `L` mean change-points.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution, including:
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
mich_matrix <- function(
  y, lambda, Q,
  fit_intercept, fit_scale,
  L, L_auto, L_max, prior,
  tol, verbose, max_iter,
  merge_prob, restart, n_search, increment, reverse,
  log_omega_l, log_pi_l
) {
  #### set up ####
  # calculate dimensions of y
  T <- nrow(y)
  d <- ncol(y)
  T1 <- T:1

  # min prob to keep component when restarting
  keep_level <- 0.9

  # log eigenvalues
  log_lambda <- log(lambda)

  # eigenvalue parameters
  omega_l <- exp(log_omega_l)
  if (prior == "var.weighted") {
    omega <- outer(T1, omega_l * lambda)
    eigen_vals <- outer(T1, lambda) + omega
    trace_omega_omega_bar_inv <- rep(d * omega_l / (1 + omega_l), T)
    log_det_omega_omega_bar <- rep(d * (log_omega_l - log(1+omega_l)), T)
  } else if (prior == 'weighted') {
    omega <- omega_l * matrix(T1, nrow = T, ncol = d, byrow = FALSE)
    eigen_vals <- outer(T1, lambda) + omega
    trace_omega_omega_bar_inv <- rep(sum(omega_l / (lambda + omega_l)), T)
    log_det_omega_omega_bar <- rep(d * log_omega_l - sum(log(lambda + omega_l)), T)
  } else if (prior == 'independent') {
    omega <- matrix(omega_l, nrow = T, ncol = d)
    eigen_vals <- outer(T1, lambda) + omega
    trace_omega_omega_bar_inv <- omega_l * rowSums(1 / eigen_vals)
    log_det_omega_omega_bar <- d * log_omega_l - rowSums(log(eigen_vals))
  }
  mean_weights <- matrix(lambda, nrow = T, ncol = d, byrow = TRUE) / eigen_vals
  sandwich_weights <- 0.5 * mean_weights * matrix(lambda, nrow = T, ncol = d, byrow = TRUE)

  # initializing posterior parameters
  post_params <- param_init(L, T, d)

  # initialize mu_0
  if (fit_intercept) {
    mu_0 <- colMeans(y[1:ceiling(log(T)),,drop=FALSE])
  } else mu_0 <- rep(0.0, d)

  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  merge_dist <- ceiling(log(T) / 2)
  while (!merged) {
    fit <- mich_matrix_cpp(
      y, mu_0, lambda, log_lambda, Q,
      trace_omega_omega_bar_inv, log_det_omega_omega_bar,
      mean_weights, sandwich_weights,
      fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
      omega, log_pi_l,
      post_params
    )

    if (merge_dist == 0 | L <= 1) break

    fit <- multi_component_merge(
      fit, lambda, mean_weights, merge_prob, merge_dist, drop = FALSE
    )

    merge_dist <- merge_dist - 1
    post_params <- fit$post_params
    merged <- fit$merged

    if (any(fit$merges)) {
      cp <- sapply(1:L, function(l) which.max(fit$post_params[[l]]$pi_bar))
      cp <- sort(unique(c(1, cp, T+1)))
      cand_scp <- list()
      dex <- 1
      for (l in 2:length(cp)) {
        QTr <- fit$QTr[cp[l-1]:(cp[l]-1),]
        Tr <- nrow(QTr)
        QTr <- QTr - matrix(colMeans(QTr[1:(Tr %/% 2),]), nrow = Tr, ncol = d, byrow = TRUE)
        null_pr <- 1 - sqrt(1 / Tr)
        log_prob_scp <- c(rep(log(1-null_pr) - log(Tr), Tr), log(null_pr))

        fit_scp <- multi_mean_scp(
          QTr, lambda, mean_weights[(T-Tr+1):T,], sandwich_weights[(T-Tr+1):T,],
          log_prob_scp + 0.5 * c(log_det_omega_omega_bar[(T-Tr+1):T], 0)
        )

        elbo_scp <- multi_elbo_fn(
          QTr, list(fit_scp), lambda, log_lambda,
          trace_omega_omega_bar_inv[(T-Tr+1):T], log_det_omega_omega_bar[(T-Tr+1):T],
          outer(Tr:1, omega_l * lambda), matrix(log_prob_scp, ncol = 1)
        )

        if (which.max(fit_scp$pi_bar) <= Tr) {
          cand_scp[[dex]] <- list(elbo = elbo_scp, post_param = fit_scp, T = Tr, lwr = cp[l-1], upr = cp[l])
          dex <- dex + 1
        }
      }

      while (length(cand_scp) > 0 & any(fit$merges)) {
        scp_dex <- which.max(sapply(seq_along(cand_scp), function(l) cand_scp[[l]]$elbo))
        fit_dex <- which.max(fit$merges)
        fit$merges[fit_dex] <- FALSE
        post_param <- cand_scp[[scp_dex]]$post_param
        Tr <-  cand_scp[[scp_dex]]$T
        lwr <-  cand_scp[[scp_dex]]$lwr
        upr <-  cand_scp[[scp_dex]]$upr

        post_param$pi_bar <- c(rep(0, lwr-1), post_param$pi_bar, rep(0, T-upr+1))
        post_param$log_pi_bar <- c(rep(0, lwr-1), post_param$log_pi_bar, rep(0, T-upr+1))
        post_param$muTLmu_bar <- c(rep(0, lwr-1), post_param$muTLmu_bar, rep(0, T-upr+1))
        post_param$QTb_bar <- rbind(
          matrix(0, nrow = lwr-1, ncol = d),
          post_param$QTb_bar,
          matrix(0, nrow = T-upr+1, ncol = d)
        )
        post_param$QTmu_bar <- rbind(
          matrix(0, nrow = lwr-1, ncol = d),
          post_param$QTmu_bar,
          matrix(post_param$QTmu_bar[Tr,], nrow = T-upr+1, ncol = d, byrow = TRUE)
        )
        post_params[[fit_dex]] <- post_param
        cand_scp[scp_dex] <- NULL
      }
    }
    # prob_order <- order(
    #   sapply(1:L, function(l) post_params[[l]][["pi_bar"]][T+1]),
    #   decreasing = TRUE
    # )
    # post_params <- post_params[prob_order]
    # log_pi_l <- log_pi_l[, prob_order, drop = FALSE]
  }

  #### auto procedure ####
  fit <- multi_auto_search(
    y, fit, log_pi_l,
    omega, trace_omega_omega_bar_inv, log_det_omega_omega_bar,
    mean_weights, sandwich_weights,
    L_max, L_auto, fit_intercept, fit_scale,
    max_iter, tol, verbose, merge_prob,
    restart, n_search, increment
  )

  #### reverse model if reverse == TRUE ####
  if (reverse) {
    if(verbose) print("reversing model")
    y <- y[T1,]
    L <- fit$L

    # reversed residuals, variance, and intercepts
    if (L_auto) log_pi_l <- matrix(log_pi_l[,1], ncol = max(1,L), nrow = T+1)
    log_pi_l[1:T,] <- log_pi_l[T1,]
    post_params <- fit$post_params
    QTr_bar_rev <- fit$QTr
    mu_0 <- colMeans(y[1:ceiling(log(T)),,drop=FALSE])

    # reversing mean components
    for (l in seq_len(L)) {
      QTr_bar_rev_l <- QTr_bar_rev + post_params[[l]][["QTmu_bar"]] - matrix(post_params[[l]][["QTmu_bar"]][T,],  nrow = T, ncol = d, byrow = TRUE)
      post_params[[l]] <- multi_mean_scp(
        QTr_bar_rev_l[T1,], lambda, mean_weights, sandwich_weights, log_pi_l[,l] + 0.5 * c(log_det_omega_omega_bar, 0)
      )
    }

    # refit model and continue search
    fit <- mich_matrix_cpp(
      y, mu_0,
      lambda, log_lambda, Q,
      trace_omega_omega_bar_inv, log_det_omega_omega_bar,
      mean_weights, sandwich_weights,
      fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
      omega, log_pi_l,
      post_params
    )

    fit <- multi_auto_search(
      y, fit, log_pi_l,
      omega, trace_omega_omega_bar_inv, log_det_omega_omega_bar,
      mean_weights, sandwich_weights,
      L_max, L_auto, fit_intercept, fit_scale,
      max_iter, tol, verbose, merge_prob,
      restart, n_search, increment
    )
  }

  #### return model ####
  # calculate mean and post prob matrix
  fit$mu_bar <- matrix(fit$mu_0, nrow = T, ncol = d, byrow = TRUE)
  QT <- t(fit$Q)
  for (l in seq_len(fit$L)) {
    fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["QTb_bar"]] %*% QT
    fit$post_params[[l]][["QTb_bar"]] <- NULL
    fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["QTmu_bar"]] %*% QT
    fit$post_params[[l]][["QTmu_bar"]] <- NULL
    fit$mu_bar <- fit$mu_bar + fit$post_params[[l]][["mu_bar"]]
  }
  if (fit$L > 0) fit$pi_bar <- sapply(1:fit$L, function(l) fit$post_params[[l]][["pi_bar"]])

  fit$y <- y
  fit$Sigma <- fit$Q %*% diag(1/fit$lambda) %*% QT
  class(fit) <- "mich.fit"

  return(fit)
}
