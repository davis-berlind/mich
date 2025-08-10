#' Univariate Multiple Independent Change-Point (MICH) Model
#'
#' @description
#' Fits the univariate version of the MICH model for mean and variance
#' change-points. Number of change-points can either be fixed or `mich_vector()`
#' will search for the number of changes that maximizes the ELBO when
#' `(L_auto | K_auto | J_auto) == TRUE`.
#'
#' @param y A numeric vector. Length \eqn{T} vector of observations.
#' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
#'   intercept is estimated, otherwise it is assumed  that
#'   \eqn{\mu_0 = 0}.
#' @param fit_scale A logical. If `fit_scale == TRUE`, then the initial
#'   precision is estimated, otherwise it is assumed  that \eqn{\lambda_0 = 1}.
#' @param standardize A logical. If `standardize == TRUE`, then `y` is centered
#'    and rescaled before fitting.
#' @param L,K,J Integers. Respective number of mean-scp, var-scp, and
#'   meanvar-scp components included in model. If `L_auto == TRUE`,
#'   `K_auto == TRUE`, or `J_auto == TRUE` then `L`, `K`, and `J` lower bound
#'   the number of each kind of change-point in the model.
#' @param L_auto,K_auto,J_auto Logicals. If `L_auto == TRUE`,
#'   `K_auto == TRUE`, and/or `J_auto == TRUE`, then `mich_vector()` returns the
#'   \eqn{L} between `L` and `L_max`, the \eqn{K} between `K` and `K_max`,
#'   and/or the \eqn{J} between `J` and `J_max` that maximize the ELBO (see
#'   Appendix C.4 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param L_max,K_max,J_max Integers.  If `L_auto == TRUE`,
#'   `K_auto == TRUE`, or `J_auto == TRUE` then `L_max`, `K_max`, and `J_max`
#'   upper bound the number of each kind of change-point in the model.
#' @param pi_l_weighted,pi_k_weighted,pi_j_weighted Logicals. If `TRUE`, then
#'   the weighted priors specified in Appendix C.2 of Berlind, Cappello, Madrid
#'   Padilla (2025) are used.
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#' @param verbose A logical. If `verbose == TRUE` and `L_auto == FALSE`,
#'   `K_auto == FALSE`, and `J_auto == FALSE` then the value of the ELBO is
#'   printed every 5000th iteration. If `verbose == TRUE` and any of `L_auto`,
#'   `K_auto`, or `J_auto` are `TRUE`, the the value of the ELBO is printed for
#'   each \eqn{(L,K,J)} as `mich_vector()` searches for the parameterization
#'   that maximized the EBLO.
#' @param max_iter An integer. Maximum number of iterations. If ELBO does not
#'   converge before `max_iter` is reached, then `converged == FALSE` in the
#'   returned fit object.
#' @param reverse A logical. If `reverse == TRUE` then MICH is fit to
#'   \eqn{\mathbf{y}_{T:1}} and the model parameters are reversed in
#'    post-processing.
#' @param detect A scalar. The detection criteria. The \eqn{i^{\text{th}}}
#'   component of the model detects a change-point only if the posterior
#'   credible set for that component contains fewer than `detect` indices.
#' @param merge_level A scalar. A value between (0,1) for the significance level
#'   to construct credible sets at when merging. A model component is only
#'   considered to be a candidate for merging if its `merge_level`-level
#'   credible set contains fewer than `detect` indices.
#' @param merge_prob A scalar. A value between (0,1) indicating the merge
#'   criterion. If the posterior probability that two components identify the
#'   same change is greater than `merge_level`, then those components are merged
#'   (see Appendix C.3 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param restart A logical. If `restart == TRUE` and `L_auto`, `K_auto`, or
#'   `J_auto` are `TRUE`, then after `n_search` increments of `L`, `K`, and/or
#'   `J`, if the ELBO has not increased, `mich_vector()` will restart by setting
#'   the `L`, `K`, and `J` components to the null model initialization (except
#'   for the components with maximum posterior probabilities > 0.9) then refit
#'   and begin the search again.
#' @param n_search An integer. Grid search parameter. Number of times to
#'   increment `L`, `K`, and/or `J` before terminating automatic procedure when
#'   `L_auto`, `K_auto`, or `J_auto` are `TRUE`.
#' @param increment An integer. Number of components to increment `L`, `K` and
#'   `J` by when `L_auto`, `K_auto`, or `J_auto` are `TRUE`.
#' @param omega_j A scalar. Prior precision parameter for meanvar-scp components
#'   of the model.
#' @param u_j A scalar. Prior shape parameter for meanvar-scp components
#'   of the model.
#' @param v_j A scalar. Prior rate parameter for meanvar-scp components
#'   of the model.
#' @param log_pi_j  A numeric matrix. A \eqn{T \times J} matrix of prior log
#'   change-point location probabilities for each of the \eqn{J} mean-variance
#'   change-points.
#' @param omega_l A scalar. Prior precision parameter for mean-scp components of
#'   the model.
#' @param log_pi_l A numeric matrix. A \eqn{T \times L} matrix of prior log
#'   change-point location probabilities for each of the \eqn{L} mean
#'   change-points.
#' @param u_k A scalar. Prior shape parameter for var-scp components
#'   of the model.
#' @param v_k A scalar. Prior rate parameter for var-scp components
#'   of the model.
#' @param log_pi_k  A numeric matrix. A \eqn{T \times K} matrix of prior log
#'   change-point location probabilities for each of the \eqn{K} variance
#'   change-points.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution, including:
#'   * `y`: A numeric vector. Original data.
#'   * `L`,`K`,`J`: Integers. Number of mean-scp, var-scp, and meanvar-scp
#'     components included in the model.
#'   * `pi_bar`: A numeric matrix. A  \eqn{T \times L} matrix of posterior
#'     change-point location probabilites.
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
#'   * `post_params`: A list. List of posterior parameters for each mean-scp
#'     component.
#'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
#'   * `converged`: A logical. Indicates whether relative increase in the ELBO
#'     is less than `tol`.
#'
mich_vector <- function(y, fit_intercept, fit_scale, standardize,
                        J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                        pi_j_weighted, pi_l_weighted, pi_k_weighted,
                        tol, verbose, max_iter, reverse,
                        detect, merge_level, merge_prob,
                        restart, n_search, increment,
                        omega_j, u_j, v_j, log_pi_j,
                        omega_l, log_pi_l,
                        u_k, v_k, log_pi_k) {

  #### set up ####
  # Flag for new model
  refit <- FALSE

  # calculate dimensions of y
  T <- length(y)

  # min prob to keep component when restarting
  keep_level <- 0.9

  # max times to merge
  merge_counter = log(T) %/% 2

  #### standardize data ####
  if(standardize) {
    center <- stats::median(y)
    scale <- stats::IQR(y)
    y <- (y - center) / scale
  }

  #### initialize mu_0 ####
  if (fit_intercept) {
    mu_0 <- mean(y[1:ceiling(2*log(T))])
  } else mu_0 <- 0.0

  #### initialize lambda_0 ####
  if (fit_scale) {
    lambda_0 <- 1 / IQR(y[1:ceiling(2*log(T))])^2
  } else lambda_0 <- 1.0

  #### initializing posterior parameters ####

  # mean components
  pi_bar_l <- matrix(1/T, nrow = T, ncol = L)
  log_pi_bar_l <- matrix(0, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  omega_bar_l <- matrix(1.0, nrow = T, ncol = L)

  # variance components
  pi_bar_k <- matrix(1/T, nrow = T, ncol = K)
  log_pi_bar_k <- matrix(0, nrow = T, ncol = K)
  u_bar_k <- u_k + (T-1:T+1) / 2
  lgamma_u_bar_k <- lgamma(u_bar_k)
  digamma_u_bar_k <- digamma(u_bar_k)
  v_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T-1:T+1) / 2

  # mean-variance components
  pi_bar_j <- matrix(1/T, nrow = T, ncol = J)
  log_pi_bar_j <- matrix(0.0, nrow = T, ncol = J)
  b_bar_j <- matrix(0.0, nrow = T, ncol = J)
  omega_bar_j <- matrix(1.0, nrow = T, ncol = J)
  u_bar_j <- u_j + (T-1:T+1) / 2
  lgamma_u_bar_j <- lgamma(u_bar_j)
  digamma_u_bar_j <- digamma(u_bar_j)

  # initialize mean components if fitting meanvar
  if (((J_auto & !L_auto & !K_auto) | J > 0) & (L + K == 0)) {
    if (verbose) print("Initializing mean components.")
    if (J_auto) log_pi_l <- sapply(1:max(1,J), function(i) log_pi_l[,1])
    else log_pi_l <- log_pi_j

    fit <- mich_vector(y, fit_intercept, fit_scale, standardize,
                       J = 0, L = J, K,
                       J_auto = L_auto, L_auto = J_auto, K_auto,
                       J_max = 0, L_max = J_max, K_max,
                       pi_j_weighted, pi_l_weighted, pi_k_weighted,
                       tol, verbose, max_iter, reverse = FALSE,
                       detect, merge_level, merge_prob,
                       restart = FALSE, n_search, increment,
                       omega_j, u_j, v_j, log_pi_j,
                       omega_l, log_pi_l,
                       u_k, v_k, log_pi_k)

    J <- fit$L
    if (J > 0) {
      refit <- TRUE
      keep <- seq_len(J)
      if (J_auto) {
        # only keep detected changes
        cred_sets <- apply(fit$mean_model$pi_bar, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        J = sum(keep)
        log_pi_j <- sapply(1:max(1,J), function(i) log_pi_j[,1])
      }

      # store initialization
      pi_bar_j <- fit$mean_model$pi_bar[, keep, drop = FALSE]
      log_pi_bar_j <- log(pi_bar_j)
      b_bar_j <- fit$mean_model$b_bar[, keep, drop = FALSE]
      omega_bar_j <- fit$mean_model$omega_bar[, keep, drop = FALSE]

      # order by longest blocks first
      chp <- apply(pi_bar_j, 2, which.max)
      chp_order <- order(chp)
      chp <- chp[chp_order]

      pi_bar_j <- pi_bar_j[,chp_order, drop = FALSE]
      log_pi_bar_j <- log_pi_bar_j[,chp_order, drop = FALSE]
      b_bar_j <- b_bar_j[,chp_order, drop = FALSE]
      omega_bar_j <- omega_bar_j[,chp_order, drop = FALSE]

      diff_order <- order(diff(c(chp, T)))

      pi_bar_j <- pi_bar_j[,diff_order, drop = FALSE]
      log_pi_bar_j <- log_pi_bar_j[,diff_order, drop = FALSE]
      b_bar_j <- b_bar_j[,diff_order, drop = FALSE]
      omega_bar_j <- omega_bar_j[,diff_order, drop = FALSE]
    }
  }
  v_bar_j <- sapply(1:max(1,J), function(i) u_bar_j)

  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  while (!merged) {
    fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                    fit_intercept, fit_scale, refit,
                    max_iter, verbose = (verbose & sum(J_auto, L_auto, K_auto) == 0), tol,
                    omega_j, u_j, v_j, log_pi_j,
                    pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j,
                    u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                    omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l,
                    u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                    lgamma_u_bar_k, digamma_u_bar_k)

    merged <- TRUE
    refit <- TRUE

    # identify components to merge
    if (L > 1) {
      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE

      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0

      mu_bar_l <- fit$mean_model$mu_bar
      mu2_bar_l <- fit$mean_model$mu2_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta

      # merge mean components with pairwise merge probabilities > merge_prob
      while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        L <- L - 1

        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        # remove model parameters set to be merged
        mu_bar_merge <- rowSums(mu_bar_l[,merge_dex])
        merge_residual <- merge_residual + mu_bar_merge
        merge_delta <- merge_delta + rowSums(mu_bar_l[,merge_dex]^2 - mu2_bar_l[,merge_dex])

        # fit mean-scp to merged residual
        merge_fit <- mean_scp(merge_residual, merge_lambda, omega_l, log_pi_l[,merge_dex[1]])

        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
          pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
          log_pi_bar_l[,merge_dex[1]] <- log_pi_bar_l[,merge_dex[2]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }

        # store merged parameters
        b_bar_l[,merge_dex[1]] <- b_bar_l[,merge_dex[1]] + b_bar_l[,merge_dex[2]]
        omega_bar_l[,merge_dex[1]] <- merge_fit$omega_bar

        # calculate merged mean signals
        mu_bar_l[,merge_dex[1]] <- mu_bar_fn(b_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])
        mu2_bar_l[,merge_dex[1]] <- mu2_bar_fn(b_bar_l[,merge_dex[1]], omega_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])

        # calculate merged residual
        merge_residual <- merge_residual - mu_bar_l[,merge_dex[1]]
        merge_delta <- merge_delta - mu_bar_l[,merge_dex[1]]^2 + mu2_bar_l[,merge_dex[1]]

        # drop merged component
        merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE]
        log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
        b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
        omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
        if (L_auto) {
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }
        mu_bar_l <- mu_bar_l[,-merge_dex[2], drop=FALSE]
        mu2_bar_l <- mu2_bar_l[,-merge_dex[2], drop=FALSE]
      }
    }

    if (K > 1) {
      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = K, nrow = K)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE

      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
      diag(merge_prob_mat) <- 0

      lambda_bar_k <- fit$var_model$lambda_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta

      # merge var components with pairwise merge probabilities > merge_prob
      while (K > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        K <- K - 1

        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        # remove model parameters set to be merged
        merge_lambda <- merge_lambda / apply(lambda_bar_k[,merge_dex], 1, prod)

        v_merge <- v_k + revcumsum(0.5 * merge_lambda * merge_delta)
        log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

        # fit var-scp to merged residual
        merge_fit <- var_scp(merge_residual, merge_lambda, u_bar_k, lgamma_u_bar_k,
                             v_merge, log_pi_merge - max(log_pi_merge))

        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_k[,merge_dex[2]]) > max(pi_bar_k[,merge_dex[1]])) {
          pi_bar_k[,merge_dex[1]] <- pi_bar_k[,merge_dex[2]]
          log_pi_bar_k[,merge_dex[1]] <- log_pi_bar_k[,merge_dex[2]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }

        # store merged parameters
        v_bar_k[,merge_dex[1]] <- merge_fit$v_bar

        # calculate merged var signal
        lambda_bar_k[,merge_dex[1]] <- lambda_bar_fn(u_bar_k, v_bar_k[,merge_dex[1]],  pi_bar_k[,merge_dex[1]])

        # calculate merged residual
        merge_lambda <- merge_lambda * lambda_bar_k[,merge_dex[1]]

        # drop merged component
        merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        pi_bar_k <- pi_bar_k[,-merge_dex[2], drop=FALSE]
        log_pi_bar_k <- log_pi_bar_k[,-merge_dex[2], drop=FALSE]
        v_bar_k <- v_bar_k[, -merge_dex[2], drop=FALSE]
        if (K_auto) {
          log_pi_k <- log_pi_k[,-merge_dex[2], drop=FALSE]
        }
        lambda_bar_k <- lambda_bar_k[,-merge_dex[2], drop=FALSE]
      }
    }

    if (J > 1) {
      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = J, nrow = J)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE

      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
      diag(merge_prob_mat) <- 0

      mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
      mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
      lambda_bar_j <- fit$meanvar_model$lambda_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta

      # merge meanvar components with pairwise merge probabilities > merge_prob
      while (J > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        J <- J - 1

        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        # remove model parameters set to be merged
        mu_bar_merge <- rowSums(mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
        merge_residual <- merge_residual + mu_bar_merge
        merge_lambda <- merge_lambda / apply(lambda_bar_j[,merge_dex], 1, prod)
        merge_delta <- merge_delta + rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
        merge_delta <- merge_delta - rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
        v_merge <- v_j + revcumsum(0.5 * merge_lambda * merge_delta)
        log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

        # fit meanvar-scp to merged residual
        merge_fit <- meanvar_scp(merge_residual, merge_lambda, omega_j, u_bar_j,
                                 lgamma_u_bar_j, v_merge, log_pi_merge - max(log_pi_merge))

        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_j[,merge_dex[2]]) > max(pi_bar_j[,merge_dex[1]])) {
          pi_bar_j[,merge_dex[1]] <- pi_bar_j[,merge_dex[2]]
          log_pi_bar_j[,merge_dex[1]] <- log_pi_bar_j[,merge_dex[2]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }

        # store merged parameters
        b_bar_j[,merge_dex[1]] <- b_bar_j[,merge_dex[1]] + b_bar_j[,merge_dex[2]]
        omega_bar_j[,merge_dex[1]] <- merge_fit$omega_bar
        v_bar_j[,merge_dex[1]] <- merge_fit$v_bar

        mu_lambda_bar_j[,merge_dex[1]] <- mu_lambda_fn(b_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
        mu2_lambda_bar_j[,merge_dex[1]] <- mu2_lambda_fn(b_bar_j[,merge_dex[1]], omega_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
        lambda_bar_j[,merge_dex[1]] <- lambda_bar_fn(u_bar_j,  v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])

        # calculate merged meanvar signals
        mu_bar_merge <- mu_lambda_bar_j[,merge_dex[1]] / lambda_bar_j[,merge_dex[1]]
        merge_lambda <- merge_lambda * lambda_bar_j[,merge_dex[1]]

        # calculate merged residual
        merge_residual <- merge_residual - mu_bar_merge
        merge_delta <- merge_delta - rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
        merge_delta <- merge_delta + rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])

        # drop merged component
        merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        pi_bar_j <- pi_bar_j[,-merge_dex[2], drop=FALSE]
        log_pi_bar_j <- log_pi_bar_j[,-merge_dex[2], drop=FALSE]
        b_bar_j <- b_bar_j[,-merge_dex[2], drop=FALSE]
        omega_bar_j <- omega_bar_j[,-merge_dex[2], drop=FALSE]
        v_bar_j <- v_bar_j[, -merge_dex[2], drop=FALSE]
        if (J_auto) {
          log_pi_j <- log_pi_j[,-merge_dex[2], drop=FALSE]
        }
        mu_lambda_bar_j <- mu_lambda_bar_j[,-merge_dex[2], drop=FALSE]
        mu2_lambda_bar_j <- mu2_lambda_bar_j[,-merge_dex[2], drop=FALSE]
        lambda_bar_j <- lambda_bar_j[,-merge_dex[2], drop=FALSE]
      }
    }
    if (verbose & !merged) print(paste0("Merging to (L = ", L, ", K = ", K, ", J = ", J, ")"))
  }

  # if components were merged out use auto procedure to increase to desired LKJ
  merge_flag <- (J < J_max & !J_auto) | (L < L_max & !L_auto) | (K < K_max & !K_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components

  #### auto procedures ####

  #### auto procedure with single component ####
  # 1. Increase L, K or J, if ELBO increases, set restart == TRUE and increase L/K/J again
  # 2. If ELBO decreases, decrease counter and fit another component
  # 3. If ELBO decreases, counter == 0, and restart == TRUE, reset components
  #    to null model (except for this with concentrated probs) and refit
  # 4. If ELBO decreases, counter == 0, restart == FALSE, and merge_counter > 0,
  #    set fit to best model so far and merge components if no merges return
  #    model, otherwise decrease merge_count and go back to 1
  # 5. If merge_count = 0 return fit

  last_restart <- ifelse(restart, 2, Inf)

  if (sum(J_auto, L_auto, K_auto) == 1 | merge_flag) {

    refit <- (L > 0 | K > 0 | J > 0)
    counter <- n_search # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo

    if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ",
                              elbo_new, "; Counter: ", counter))

    # continue search until n_search exhausted or max components exceeded
    while (J < J_max | L < L_max | K < K_max) {
      if (J_auto | J < J_max) {
        # increment dimension of parameters
        J <- J + increment
        if (J > 1 & J_auto) {
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }
        pi_bar_j <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_j)
        log_pi_bar_j <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_j)
        b_bar_j <- cbind(matrix(0.0, nrow = T, ncol = increment), b_bar_j)
        omega_bar_j <- cbind(matrix(1.0, nrow = T, ncol = increment), omega_bar_j)
        v_bar_j <- cbind(matrix(v_j + (T-1:T+1) / 2, nrow = T, ncol = increment), v_bar_j)
      }

      if (L_auto | L < L_max) {
        # increment dimension of parameters
        L <- L + increment
        if (L > 1 & L_auto) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }
        pi_bar_l <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_l)
        log_pi_bar_l <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_l)
        b_bar_l <- cbind(matrix(0.0, nrow = T, ncol = increment), b_bar_l)
        omega_bar_l <- cbind(matrix(1.0, nrow = T, ncol = increment), omega_bar_l)
      }
      if (K_auto | K < K_max) {
        # increment dimension of parameters
        K <- K + increment
        if (K > 1 & K_auto) {
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }
        pi_bar_k <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_k)
        log_pi_bar_k <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_k)
        v_bar_k <- cbind(matrix(v_k + (T-1:T+1) / 2, nrow = T, ncol = increment), v_bar_k)
      }

      # fit incremented model
      fit_new <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                          fit_intercept, fit_scale, refit,
                          max_iter = max_iter,
                          verbose = FALSE, tol,
                          omega_j, u_j, v_j, log_pi_j,
                          pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j,
                          u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                          omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l,
                          u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                          lgamma_u_bar_k, digamma_u_bar_k)

      # test if model improved or merge/restart
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new, "; Counter: ", counter))

      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE
      } else if (counter == 0 & restart) {
        if (verbose) print(paste0("Restarting at (L = ", L, ", K = ", K, ", J = ", J, ")"))
        restart <- FALSE
        counter <- n_search

        if (J_auto & J > last_restart) {
          last_restart <- J

          # reorder by longest blocks first
          chp <- apply(pi_bar_j, 2, which.max)
          chp_order <- order(chp)
          chp <- chp[chp_order]

          pi_bar_j <- pi_bar_j[,chp_order, drop = FALSE]
          log_pi_bar_j <- log_pi_bar_j[,chp_order, drop = FALSE]
          b_bar_j <- b_bar_j[,chp_order, drop = FALSE]
          omega_bar_j <- omega_bar_j[,chp_order, drop = FALSE]
          v_bar_j <- v_bar_j[,chp_order, drop = FALSE]

          diff_order <- order(diff(c(chp, T)))

          pi_bar_j <- pi_bar_j[,diff_order, drop = FALSE]
          log_pi_bar_j <- log_pi_bar_j[,diff_order, drop = FALSE]
          b_bar_j <- b_bar_j[,diff_order, drop = FALSE]
          omega_bar_j <- omega_bar_j[,diff_order, drop = FALSE]
          v_bar_j <- v_bar_j[,diff_order, drop = FALSE]

          # identify components with max prob > keep_level
          keep <- apply(pi_bar_j, 2, max) > keep_level

          # reset components with max prob < keep_level to null components
          if (sum(keep) < J) {
            keep_inc <- max(sum(!keep) - increment, 0)
            J <- sum(keep) + keep_inc

            log_pi_j <- sapply(1:J, function(i) log_pi_j[,1, drop = FALSE])
            pi_bar_j <- cbind(matrix(1/T, nrow = T, ncol = keep_inc), pi_bar_j[,keep, drop = FALSE])
            log_pi_bar_j <- cbind(matrix(0.0, nrow = T, ncol = keep_inc), log_pi_bar_j[,keep, drop = FALSE])
            b_bar_j <- cbind(matrix(0.0, nrow = T, ncol = keep_inc), b_bar_j[,keep, drop = FALSE])
            omega_bar_j <- cbind(matrix(1.0, nrow = T, ncol = keep_inc), omega_bar_j[,keep, drop = FALSE])
            v_bar_j <- cbind(matrix(u_j, nrow = T, ncol = keep_inc, byrow = TRUE) + (T-1:T+1) / 2, v_bar_j[,keep, drop = FALSE])
          }
        }

        if (L_auto & L > last_restart) {
          last_restart <- L

          # reorder by longest blocks first
          chp <- apply(pi_bar_l, 2, which.max)
          chp_order <- order(chp)
          chp <- chp[chp_order]

          pi_bar_l <- pi_bar_l[,chp_order, drop = FALSE]
          log_pi_bar_l <- log_pi_bar_l[,chp_order, drop = FALSE]
          b_bar_l <- b_bar_l[,chp_order, drop = FALSE]
          omega_bar_l <- omega_bar_l[,chp_order, drop = FALSE]

          diff_order <- order(diff(c(chp, T)))

          pi_bar_l <- pi_bar_l[,diff_order, drop = FALSE]
          log_pi_bar_l <- log_pi_bar_l[,diff_order, drop = FALSE]
          b_bar_l <- b_bar_l[,diff_order, drop = FALSE]
          omega_bar_l <- omega_bar_l[,diff_order, drop = FALSE]

          # identify components with max prob > keep_level
          keep <- apply(pi_bar_l, 2, max) > keep_level

          # reset components with max prob < keep_level to null components
          if (sum(keep) < L) {
            keep_inc <- max(sum(!keep) - increment, 0)
            L <- sum(keep) + keep_inc

            log_pi_l <- sapply(1:L, function(i) log_pi_l[,1, drop = FALSE])
            pi_bar_l <- cbind(matrix(1/T, nrow = T, ncol =keep_inc), pi_bar_l[,keep, drop = FALSE])
            log_pi_bar_l <- cbind(matrix(0.0, nrow = T, ncol = keep_inc), log_pi_bar_l[,keep, drop = FALSE])
            b_bar_l <- cbind(matrix(0.0, nrow = T, ncol = keep_inc), b_bar_l[,keep, drop = FALSE])
            omega_bar_l <- cbind(matrix(1.0, nrow = T, ncol = keep_inc), omega_bar_l[,keep, drop = FALSE])
          }
        }

        if (K_auto & K > last_restart) {
          last_restart <- K

          # reorder by longest blocks first
          chp <- apply(pi_bar_k, 2, which.max)
          chp_order <- order(chp)
          chp <- chp[chp_order]

          pi_bar_k <- pi_bar_k[,chp_order]
          log_pi_bar_k <- log_pi_bar_k[,chp_order]
          v_bar_k <- v_bar_k[,chp_order]

          diff_order <- order(diff(c(chp, T)))

          pi_bar_k <- pi_bar_k[,diff_order]
          log_pi_bar_k <- log_pi_bar_k[,diff_order]
          v_bar_k <- v_bar_k[,diff_order]

          # identify components with max prob > keep_level
          keep <- apply(pi_bar_k, 2, max) > keep_level

          # reset components with max prob < keep_level to null components

                    if (sum(keep) < K) {
            keep_inc <- max(sum(!keep) - increment, 0)
            K <- sum(keep) + keep_inc

            log_pi_k <- sapply(1:K, function(i) log_pi_k[,1, drop = FALSE])
            pi_bar_k <- cbind(matrix(1/T, nrow = T, ncol = keep_inc), pi_bar_k[,keep, drop = FALSE])
            log_pi_bar_k <- cbind(matrix(0.0, nrow = T, ncol = keep_inc), log_pi_bar_k[,keep, drop = FALSE])
            v_bar_k <- cbind(matrix(u_k, nrow = T, ncol = keep_inc, byrow = TRUE) + (T-1:T+1) / 2, v_bar_k[,keep, drop = FALSE])
          }
        }
      } else if (counter == 0 | (L == L_max & K == K_max & J == J_max)) {
        # merging components
        counter <- n_search # reset counter in case searching continues
        merge_counter <- merge_counter - 1
        no_merges <- TRUE  # flag indicating no components were merged
        merged <- FALSE # flag indicating components have been merged
        if (verbose) print(paste0("Merging. Merge Counter: ", merge_counter))

        # set fit to best model so far
        L <- fit$L; K <- fit$K; J <- fit$J;
        merge_residual <- fit$residual
        merge_lambda <- fit$lambda
        merge_delta <- fit$delta

        if (L >= 1) {
          log_pi_l <- log_pi_l[, 1:L, drop = FALSE]
          pi_bar_l <- fit$mean_model$pi_bar
          log_pi_bar_l <- log(pi_bar_l)
          b_bar_l <- fit$mean_model$b_bar
          omega_bar_l <- fit$mean_model$omega_bar
          mu_bar_l <- fit$mean_model$mu_bar
          mu2_bar_l <- fit$mean_model$mu2_bar
        }

        if (K >= 1) {
          log_pi_k <- log_pi_k[, 1:K, drop = FALSE]
          pi_bar_k <- fit$var_model$pi_bar
          log_pi_bar_k <- log(pi_bar_k)
          v_bar_k <- fit$var_model$v_bar
          lambda_bar_k <- fit$var_model$lambda_bar
        }

        if (J >= 1) {
          log_pi_j <- log_pi_j[, 1:J, drop = FALSE]
          pi_bar_j <- fit$meanvar_model$pi_bar
          log_pi_bar_j <- log(pi_bar_j)
          b_bar_j <- fit$meanvar_model$b_bar
          omega_bar_j <- fit$meanvar_model$omega_bar
          v_bar_j <- fit$meanvar_model$v_bar
          mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
          mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
          lambda_bar_j <- fit$meanvar_model$lambda_bar
        }

        while (!merged) {
          fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                          fit_intercept, fit_scale, refit = TRUE,
                          max_iter = ifelse(no_merges, 1, max_iter),
                          verbose = FALSE, tol,
                          omega_j, u_j, v_j, log_pi_j,
                          pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j,
                          u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                          omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l,
                          u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                          lgamma_u_bar_k, digamma_u_bar_k)

          merged <- TRUE

          # identify components to merge
          if (L > 1) {
            # only merge columns with credible sets with length less than detect
            cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = L, nrow = L)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE

            # compute pairwise merge probabilities
            merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
            diag(merge_prob_mat) <- 0

            mu_bar_l <- fit$mean_model$mu_bar
            mu2_bar_l <- fit$mean_model$mu2_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta

            # merge mean components with pairwise merge probabilities > merge_prob
            while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
              merged <- FALSE
              L <- L - 1

              # identify components with largest pairwise merge probabilities
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

              # remove model parameters set to be merged
              mu_bar_merge <- rowSums(mu_bar_l[,merge_dex])
              merge_residual <- merge_residual + mu_bar_merge
              merge_delta <- merge_delta + rowSums(mu_bar_l[,merge_dex]^2 - mu2_bar_l[,merge_dex])

              # fit mean-scp to merged residual
              merge_fit <- mean_scp(merge_residual, merge_lambda, omega_l, log_pi_l[,merge_dex[1]])

              # keep probabilities of component with largest posterior probability
              if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
                pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
                log_pi_bar_l[,merge_dex[1]] <- log_pi_bar_l[,merge_dex[2]]
                merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
              }

              # store merged parameters
              b_bar_l[,merge_dex[1]] <- b_bar_l[,merge_dex[1]] + b_bar_l[,merge_dex[2]]
              omega_bar_l[,merge_dex[1]] <- merge_fit$omega_bar

              # calculate merged mean signals
              mu_bar_l[,merge_dex[1]] <- mu_bar_fn(b_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])
              mu2_bar_l[,merge_dex[1]] <- mu2_bar_fn(b_bar_l[,merge_dex[1]], omega_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])

              # calculate merged residual
              merge_residual <- merge_residual - mu_bar_l[,merge_dex[1]]
              merge_delta <- merge_delta - mu_bar_l[,merge_dex[1]]^2 + mu2_bar_l[,merge_dex[1]]

              # drop merged component
              merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
              keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
              pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE]
              log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
              b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
              omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
              if (L_auto) {
                log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
              }
              mu_bar_l <- mu_bar_l[,-merge_dex[2], drop=FALSE]
              mu2_bar_l <- mu2_bar_l[,-merge_dex[2], drop=FALSE]
            }
          }

          if (K > 1) {
            # only merge columns with credible sets with length less than detect
            cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = K, nrow = K)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE

            # compute pairwise merge probabilities
            merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
            diag(merge_prob_mat) <- 0

            lambda_bar_k <- fit$var_model$lambda_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta

            # merge var components with pairwise merge probabilities > merge_prob
            while (K > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
              merged <- FALSE
              K <- K - 1

              # identify components with largest pairwise merge probabilities
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

              # remove model parameters set to be merged
              merge_lambda <- merge_lambda / apply(lambda_bar_k[,merge_dex], 1, prod)

              v_merge <- v_k + revcumsum(0.5 * merge_lambda * merge_delta)
              log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

              # fit var-scp to merged residual
              merge_fit <- var_scp(merge_residual, merge_lambda, u_bar_k, lgamma_u_bar_k,
                                   v_merge, log_pi_merge - max(log_pi_merge))

              # keep probabilities of component with largest posterior probability
              if (max(pi_bar_k[,merge_dex[2]]) > max(pi_bar_k[,merge_dex[1]])) {
                pi_bar_k[,merge_dex[1]] <- pi_bar_k[,merge_dex[2]]
                log_pi_bar_k[,merge_dex[1]] <- log_pi_bar_k[,merge_dex[2]]
                merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
              }

              # store merged parameters
              v_bar_k[,merge_dex[1]] <- merge_fit$v_bar

              # calculate merged var signal
              lambda_bar_k[,merge_dex[1]] <- lambda_bar_fn(u_bar_k, v_bar_k[,merge_dex[1]],  pi_bar_k[,merge_dex[1]])

              # calculate merged residual
              merge_lambda <- merge_lambda * lambda_bar_k[,merge_dex[1]]

              # drop merged component
              merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
              keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
              pi_bar_k <- pi_bar_k[,-merge_dex[2], drop=FALSE]
              log_pi_bar_k <- log_pi_bar_k[,-merge_dex[2], drop=FALSE]
              v_bar_k <- v_bar_k[, -merge_dex[2], drop=FALSE]
              if (K_auto) {
                log_pi_k <- log_pi_k[,-merge_dex[2], drop=FALSE]
              }
              lambda_bar_k <- lambda_bar_k[,-merge_dex[2], drop=FALSE]
            }
          }

          if (J > 1) {
            # only merge columns with credible sets with length less than detect
            cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = J, nrow = J)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE

            # compute pairwise merge probabilities
            merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
            diag(merge_prob_mat) <- 0

            mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
            mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
            lambda_bar_j <- fit$meanvar_model$lambda_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta

            # merge meanvar components with pairwise merge probabilities > merge_prob
            while (J > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
              merged <- FALSE
              J <- J - 1

              # identify components with largest pairwise merge probabilities
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

              # remove model parameters set to be merged
              mu_bar_merge <- rowSums(mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
              merge_residual <- merge_residual + mu_bar_merge
              merge_lambda <- merge_lambda / apply(lambda_bar_j[,merge_dex], 1, prod)
              merge_delta <- merge_delta + rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
              merge_delta <- merge_delta - rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
              v_merge <- v_j + revcumsum(0.5 * merge_lambda * merge_delta)
              log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

              # fit meanvar-scp to merged residual
              merge_fit <- meanvar_scp(merge_residual, merge_lambda, omega_j, u_bar_j,
                                       lgamma_u_bar_j, v_merge, log_pi_merge - max(log_pi_merge))

              # keep probabilities of component with largest posterior probability
              if (max(pi_bar_j[,merge_dex[2]]) > max(pi_bar_j[,merge_dex[1]])) {
                pi_bar_j[,merge_dex[1]] <- pi_bar_j[,merge_dex[2]]
                log_pi_bar_j[,merge_dex[1]] <- log_pi_bar_j[,merge_dex[2]]
                merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
              }

              # store merged parameters
              b_bar_j[,merge_dex[1]] <- b_bar_j[,merge_dex[1]] + b_bar_j[,merge_dex[2]]
              omega_bar_j[,merge_dex[1]] <- merge_fit$omega_bar
              v_bar_j[,merge_dex[1]] <- merge_fit$v_bar

              mu_lambda_bar_j[,merge_dex[1]] <- mu_lambda_fn(b_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
              mu2_lambda_bar_j[,merge_dex[1]] <- mu2_lambda_fn(b_bar_j[,merge_dex[1]], omega_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
              lambda_bar_j[,merge_dex[1]] <- lambda_bar_fn(u_bar_j,  v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])

              # calculate merged meanvar signals
              mu_bar_merge <- mu_lambda_bar_j[,merge_dex[1]] / lambda_bar_j[,merge_dex[1]]
              merge_lambda <- merge_lambda * lambda_bar_j[,merge_dex[1]]

              # calculate merged residual
              merge_residual <- merge_residual - mu_bar_merge
              merge_delta <- merge_delta - rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
              merge_delta <- merge_delta + rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])

              # drop merged component
              merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
              keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
              pi_bar_j <- pi_bar_j[,-merge_dex[2], drop=FALSE]
              log_pi_bar_j <- log_pi_bar_j[,-merge_dex[2], drop=FALSE]
              b_bar_j <- b_bar_j[,-merge_dex[2], drop=FALSE]
              omega_bar_j <- omega_bar_j[,-merge_dex[2], drop=FALSE]
              v_bar_j <- v_bar_j[, -merge_dex[2], drop=FALSE]
              if (J_auto) {
                log_pi_j <- log_pi_j[,-merge_dex[2], drop=FALSE]
              }
              mu_lambda_bar_j <- mu_lambda_bar_j[,-merge_dex[2], drop=FALSE]
              mu2_lambda_bar_j <- mu2_lambda_bar_j[,-merge_dex[2], drop=FALSE]
              lambda_bar_j <- lambda_bar_j[,-merge_dex[2], drop=FALSE]
            }
          }
          if (verbose & !merged) print(paste0("Merging to (L = ", L, ", K = ", K, ", J = ", J, ")"))
        }

        elbo <- max(fit$elbo)
        if (elbo > merge_elbo) {
          merge_elbo <- elbo
          fit_merge <- fit
        }

        if (no_merges | merge_counter == 0) {
          fit <- fit_merge
          break
        }
      } else {
        counter <- counter - 1
      }
    }
  }

  #### auto procedure with multiple components ####
  # 0. Set counter = (J_auto + L_auto + K_auto) * n_search
  # 1. Separately increment L, K, J and refit model, go in direction of max ELBO
  #    increase, if ELBO increases reset counter otherwise decrement counter.
  # 2. If counter == 0 return model with largest ELBO.
  if (sum(J_auto, L_auto, K_auto) > 1) {

    auto_done <- FALSE
    counter <- (J_auto + L_auto + K_auto) * n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo))

    while (TRUE) {
      elbo_new <- rep(-Inf, 3)
      if (J_auto) {
        J <- J + increment
        if (J > 1) {
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }

        # fit new model
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_search, increment,
                               omega_j, u_j, v_j, log_pi_j,
                               omega_l, log_pi_l,
                               u_k, v_k, log_pi_k)

        elbo_new[1] <- fit_new$elbo[length(fit_new$elbo)]

        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[1]))
        if (elbo_new[1] > elbo) {
          fit <- fit_new
          elbo <- elbo_new[1]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        J <- J - increment
        if (J > 0) {
          log_pi_j <- sapply(1:J, function(i) log_pi_j[,1])
        }
      }

      if (L_auto) {
        L <- L + increment
        if (L > 1) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }

        # fit new model
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_search, increment,
                               omega_j, u_j, v_j, log_pi_j,
                               omega_l, log_pi_l,
                               u_k, v_k, log_pi_k)

        elbo_new[2] <- fit_new$elbo[length(fit_new$elbo)]

        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[2]))
        if (elbo_new[2] > elbo) {
          fit <- fit_new
          elbo <- elbo_new[2]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        L <- L - increment
        if (L > 0) {
          log_pi_l <- sapply(1:L, function(i) log_pi_l[,1])
        }
      }

      if (K_auto) {
        K <- K + increment
        if (K > 1) {
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }

        # fit new model
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_search,
                               omega_j, u_j, v_j, log_pi_j,
                               omega_l, log_pi_l,
                               u_k, v_k, log_pi_k)

        elbo_new[3] <- fit_new$elbo[length(fit_new$elbo)]

        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[3]))
        if (elbo_new[3] > elbo) {
          fit <- fit_new
          elbo <- elbo_new[3]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        K <- K - 1
        if (K > 0) {
          log_pi_k <- sapply(1:K, function(i) log_pi_k[,1])
        }
      }
      if (counter <= 0) {
        auto_done <- TRUE
        break
      }

      if (which.max(elbo_new) == 1) {
        J <- J + increment
        if (J > 1) {
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }
      } else if (which.max(elbo_new) == 2) {
        L <- L + increment
        if (L > 1) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }
      } else {
        K <- K + increment
        if (K > 1) {
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }
      }
      if (auto_done) break
    }
  }

  #### reverse model if reverse == TRUE ####
  if (reverse){
    if(verbose) print("reversing model")

    L <- fit$L; K <- fit$K; J <- fit$J;

    # reversed residuals, variance, and intercepts
    r_tilde <- fit$residual[T:1]
    lambda_bar <- fit$lambda[T:1]
    delta <- fit$delta[T:1]
    mu_0 <- fit$mu[T]
    lambda_0 <- lambda_bar[1]

    # don't reverse weighted priors
    if (!pi_l_weighted) {
      log_pi_l <- log_pi_l[T:1,1:max(1,L),drop = FALSE]
    } else {
      log_pi_l <- sapply(1:max(1,L), function(i) log_pi_l[,1])
    }
    if (!pi_k_weighted) {
      log_pi_k <- log_pi_k[T:1,1:max(1,K),drop = FALSE]
    } else {
      log_pi_k <- sapply(1:max(1,K), function(i) log_pi_k[,1])
    }
    if (!pi_j_weighted) {
      log_pi_j <- log_pi_j[T:1,1:max(1,J),drop = FALSE]
    } else {
      log_pi_j <- sapply(1:max(1,J), function(i) log_pi_j[,1])
    }

    # reversing mean components
    L_seq = seq_len(L)
    pi_bar_l <- matrix(nrow = T, ncol = L)
    log_pi_bar_l <- matrix(nrow = T, ncol = L)
    b_bar_l <- matrix(nrow = T, ncol = L)
    omega_bar_l <- matrix(nrow = T, ncol = L)

    for (l in L_seq) {
      tau_l <- which.max(fit$mean_model$pi_bar[,l])
      mu_bar_l <- fit$mean_model$mu_bar[,l]
      r_tilde_l <- fit$residual + mu_bar_l - mean(mu_bar_l[tau_l:T])
      r_tilde_l <- r_tilde_l[T:1]

      fit_scp <- mean_scp(r_tilde_l, lambda_bar, omega_l, log_pi_l[,l])

      pi_bar_l[,l] <- fit_scp$pi_bar
      log_pi_bar_l[,l] <- fit_scp$log_pi_bar
      b_bar_l[,l] <- fit_scp$b_bar
      omega_bar_l[,l] <- fit_scp$omega_bar
    }

    # reversing variance components
    K_seq = seq_len(K)
    pi_bar_k <- matrix(nrow = T, ncol = K)
    log_pi_bar_k <- matrix(nrow = T, ncol = K)
    v_bar_k <- matrix(nrow = T, ncol = K)

    for (k in K_seq) {
      tau_k <- which.max(fit$var_model$pi_bar[,k])
      lambda_bar_k <- fit$lambda / (fit$var_model$lambda_bar[,k] / mean(fit$var_model$lambda_bar[tau_k:T, k]))
      lambda_bar_k <- lambda_bar_k[T:1]

      fit_scp <- var_scp(r_tilde, lambda_bar_k, u_bar_k, lgamma_u_bar_k,
                         rep(v_k, T), log_pi_k[,k])

      pi_bar_k[,k] <- fit_scp$pi_bar
      log_pi_bar_k[,k] <- fit_scp$log_pi_bar
      v_bar_k[,k] <- fit_scp$v_bar
    }

    # reversing mean-variance components
    J_seq = seq_len(J)
    pi_bar_j <- matrix(nrow = T, ncol = J)
    log_pi_bar_j <- matrix(nrow = T, ncol = J)
    b_bar_j <- matrix(nrow = T, ncol = J)
    omega_bar_j <- matrix(nrow = T, ncol = J)
    v_bar_j <- matrix(nrow = T, ncol = J)

    for (j in J_seq) {
      tau_j <- which.max(fit$meanvar_model$pi_bar[,j])
      mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar[,j]
      mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar[,j]
      lambda_bar_j <- fit$meanvar_model$lambda_bar[,j]
      r_tilde_j <- fit$residual + mu_lambda_bar_j / lambda_bar_j - mean(mu_lambda_bar_j[tau_j:T] / lambda_bar_j[tau_j:T])
      r_tilde_j <- r_tilde_j[T:1]
      lambda_bar_j <- fit$lambda / (fit$meanvar_model$lambda_bar[,j] / mean(fit$meanvar_model$lambda_bar[tau_j:T,j]))
      lambda_bar_j <- lambda_bar_j[T:1]

      fit_scp <- meanvar_scp(r_tilde_j, lambda_bar_j, omega_j, u_bar_j,
                             lgamma_u_bar_j, rep(v_j, T), log_pi_j[,j])

      pi_bar_j[,j] <- fit_scp$pi_bar
      log_pi_bar_j[,j] <- fit_scp$log_pi_bar
      b_bar_j[,j] <- fit_scp$b_bar
      omega_bar_j[,j] <- fit_scp$omega_bar
      v_bar_j[,j] <- fit_scp$v_bar
    }

    # refit model and merge components
    merged <- FALSE  # flag indicating components have been merged
    while (!merged) {
      fit <- mich_cpp(y[T:1], J, L, K, mu_0, lambda_0,
                      fit_intercept, fit_scale, refit = TRUE,
                      max_iter, verbose = FALSE, tol,
                      omega_j, u_j, v_j, log_pi_j,
                      pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j,
                      u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                      omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l,
                      u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                      lgamma_u_bar_k, digamma_u_bar_k)

      merged <- TRUE

      # identify components to merge
      if (L > 1) {
        # only merge columns with credible sets with length less than detect
        cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        keep_mat <- matrix(FALSE, ncol = L, nrow = L)
        keep_mat[keep, keep] <- TRUE
        diag(keep_mat) <- FALSE

        # compute pairwise merge probabilities
        merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
        diag(merge_prob_mat) <- 0

        mu_bar_l <- fit$mean_model$mu_bar
        mu2_bar_l <- fit$mean_model$mu2_bar
        merge_residual <- fit$residual
        merge_lambda <- fit$lambda
        merge_delta <- fit$delta

        # merge mean components with pairwise merge probabilities > merge_prob
        while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
          merged <- FALSE
          L <- L - 1

          # identify components with largest pairwise merge probabilities
          merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

          # remove model parameters set to be merged
          mu_bar_merge <- rowSums(mu_bar_l[,merge_dex])
          merge_residual <- merge_residual + mu_bar_merge
          merge_delta <- merge_delta + rowSums(mu_bar_l[,merge_dex]^2 - mu2_bar_l[,merge_dex])

          # fit mean-scp to merged residual
          merge_fit <- mean_scp(merge_residual, merge_lambda, omega_l, log_pi_l[,merge_dex[1]])

          # keep probabilities of component with largest posterior probability
          if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
            pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
            log_pi_bar_l[,merge_dex[1]] <- log_pi_bar_l[,merge_dex[2]]
            merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
            merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
          }

          # store merged parameters
          b_bar_l[,merge_dex[1]] <- b_bar_l[,merge_dex[1]] + b_bar_l[,merge_dex[2]]
          omega_bar_l[,merge_dex[1]] <- merge_fit$omega_bar

          # calculate merged mean signals
          mu_bar_l[,merge_dex[1]] <- mu_bar_fn(b_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])
          mu2_bar_l[,merge_dex[1]] <- mu2_bar_fn(b_bar_l[,merge_dex[1]], omega_bar_l[,merge_dex[1]], pi_bar_l[,merge_dex[1]])

          # calculate merged residual
          merge_residual <- merge_residual - mu_bar_l[,merge_dex[1]]
          merge_delta <- merge_delta - mu_bar_l[,merge_dex[1]]^2 + mu2_bar_l[,merge_dex[1]]

          # drop merged component
          merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
          keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
          pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE]
          log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
          b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
          omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
          if (L_auto) {
            log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
          }
          mu_bar_l <- mu_bar_l[,-merge_dex[2], drop=FALSE]
          mu2_bar_l <- mu2_bar_l[,-merge_dex[2], drop=FALSE]
        }
      }

      if (K > 1) {
        # only merge columns with credible sets with length less than detect
        cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        keep_mat <- matrix(FALSE, ncol = K, nrow = K)
        keep_mat[keep, keep] <- TRUE
        diag(keep_mat) <- FALSE

        # compute pairwise merge probabilities
        merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
        diag(merge_prob_mat) <- 0

        lambda_bar_k <- fit$var_model$lambda_bar
        merge_residual <- fit$residual
        merge_lambda <- fit$lambda
        merge_delta <- fit$delta

        # merge var components with pairwise merge probabilities > merge_prob
        while (K > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
          merged <- FALSE
          K <- K - 1

          # identify components with largest pairwise merge probabilities
          merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

          # remove model parameters set to be merged
          merge_lambda <- merge_lambda / apply(lambda_bar_k[,merge_dex], 1, prod)

          v_merge <- v_k + revcumsum(0.5 * merge_lambda * merge_delta)
          log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

          # fit var-scp to merged residual
          merge_fit <- var_scp(merge_residual, merge_lambda, u_bar_k, lgamma_u_bar_k,
                               v_merge, log_pi_merge - max(log_pi_merge))

          # keep probabilities of component with largest posterior probability
          if (max(pi_bar_k[,merge_dex[2]]) > max(pi_bar_k[,merge_dex[1]])) {
            pi_bar_k[,merge_dex[1]] <- pi_bar_k[,merge_dex[2]]
            log_pi_bar_k[,merge_dex[1]] <- log_pi_bar_k[,merge_dex[2]]
            merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
            merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
          }

          # store merged parameters
          v_bar_k[,merge_dex[1]] <- merge_fit$v_bar

          # calculate merged var signal
          lambda_bar_k[,merge_dex[1]] <- lambda_bar_fn(u_bar_k, v_bar_k[,merge_dex[1]],  pi_bar_k[,merge_dex[1]])

          # calculate merged residual
          merge_lambda <- merge_lambda * lambda_bar_k[,merge_dex[1]]

          # drop merged component
          merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
          keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
          pi_bar_k <- pi_bar_k[,-merge_dex[2], drop=FALSE]
          log_pi_bar_k <- log_pi_bar_k[,-merge_dex[2], drop=FALSE]
          v_bar_k <- v_bar_k[, -merge_dex[2], drop=FALSE]
          if (K_auto) {
            log_pi_k <- log_pi_k[,-merge_dex[2], drop=FALSE]
          }
          lambda_bar_k <- lambda_bar_k[,-merge_dex[2], drop=FALSE]
        }
      }

      if (J > 1) {
        # only merge columns with credible sets with length less than detect
        cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        keep_mat <- matrix(FALSE, ncol = J, nrow = J)
        keep_mat[keep, keep] <- TRUE
        diag(keep_mat) <- FALSE

        # compute pairwise merge probabilities
        merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
        diag(merge_prob_mat) <- 0

        mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
        mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
        lambda_bar_j <- fit$meanvar_model$lambda_bar
        merge_residual <- fit$residual
        merge_lambda <- fit$lambda
        merge_delta <- fit$delta

        # merge meanvar components with pairwise merge probabilities > merge_prob
        while (J > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
          merged <- FALSE
          J <- J - 1

          # identify components with largest pairwise merge probabilities
          merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

          # remove model parameters set to be merged
          mu_bar_merge <- rowSums(mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
          merge_residual <- merge_residual + mu_bar_merge
          merge_lambda <- merge_lambda / apply(lambda_bar_j[,merge_dex], 1, prod)
          merge_delta <- merge_delta + rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
          merge_delta <- merge_delta - rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
          v_merge <- v_j + revcumsum(0.5 * merge_lambda * merge_delta)
          log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))

          # fit meanvar-scp to merged residual
          merge_fit <- meanvar_scp(merge_residual, merge_lambda, omega_j, u_bar_j,
                                   lgamma_u_bar_j, v_merge, log_pi_merge - max(log_pi_merge))

          # keep probabilities of component with largest posterior probability
          if (max(pi_bar_j[,merge_dex[2]]) > max(pi_bar_j[,merge_dex[1]])) {
            pi_bar_j[,merge_dex[1]] <- pi_bar_j[,merge_dex[2]]
            log_pi_bar_j[,merge_dex[1]] <- log_pi_bar_j[,merge_dex[2]]
            merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
            merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
          }

          # store merged parameters
          b_bar_j[,merge_dex[1]] <- b_bar_j[,merge_dex[1]] + b_bar_j[,merge_dex[2]]
          omega_bar_j[,merge_dex[1]] <- merge_fit$omega_bar
          v_bar_j[,merge_dex[1]] <- merge_fit$v_bar

          mu_lambda_bar_j[,merge_dex[1]] <- mu_lambda_fn(b_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
          mu2_lambda_bar_j[,merge_dex[1]] <- mu2_lambda_fn(b_bar_j[,merge_dex[1]], omega_bar_j[,merge_dex[1]], u_bar_j, v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])
          lambda_bar_j[,merge_dex[1]] <- lambda_bar_fn(u_bar_j,  v_bar_j[,merge_dex[1]], pi_bar_j[,merge_dex[1]])

          # calculate merged meanvar signals
          mu_bar_merge <- mu_lambda_bar_j[,merge_dex[1]] / lambda_bar_j[,merge_dex[1]]
          merge_lambda <- merge_lambda * lambda_bar_j[,merge_dex[1]]

          # calculate merged residual
          merge_residual <- merge_residual - mu_bar_merge
          merge_delta <- merge_delta - rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2)
          merge_delta <- merge_delta + rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])

          # drop merged component
          merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
          keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
          pi_bar_j <- pi_bar_j[,-merge_dex[2], drop=FALSE]
          log_pi_bar_j <- log_pi_bar_j[,-merge_dex[2], drop=FALSE]
          b_bar_j <- b_bar_j[,-merge_dex[2], drop=FALSE]
          omega_bar_j <- omega_bar_j[,-merge_dex[2], drop=FALSE]
          v_bar_j <- v_bar_j[, -merge_dex[2], drop=FALSE]
          if (J_auto) {
            log_pi_j <- log_pi_j[,-merge_dex[2], drop=FALSE]
          }
          mu_lambda_bar_j <- mu_lambda_bar_j[,-merge_dex[2], drop=FALSE]
          mu2_lambda_bar_j <- mu2_lambda_bar_j[,-merge_dex[2], drop=FALSE]
          lambda_bar_j <- lambda_bar_j[,-merge_dex[2], drop=FALSE]
        }
      }
      if (verbose & !merged) print(paste0("Merging to (L = ", L, ", K = ", K, ", J = ", J, ")"))
    }
  }

  #### return model ####
  class(fit) <- "mich.fit"

  # rescale data to original units
  if (standardize) {
    fit$y <- fit$y * scale + center
    fit$residual <- fit$residual * scale
    fit$delta <- fit$delta * scale^2
    fit$mu_0 <- fit$mu_0 * scale + center
    fit$mu <- fit$mu * scale + center
    fit$lambda_0 <- fit$lambda_0 / scale^2
    fit$lambda <- fit$lambda / scale^2

    if (fit$L > 0) {
      fit$mean_model$b_bar <- fit$mean_model$b_bar * scale
      fit$mean_model$omega_bar <- fit$mean_model$omega_bar / scale^2
      fit$mean_model$mu_bar <- fit$mean_model$mu_bar * scale
      fit$mean_model$mu2_bar <- fit$mean_model$mu2_bar * scale^2
    }

    if (fit$J > 0) {
      fit$meanvar_model$b_bar <- fit$meanvar_model$b_bar * scale
      fit$meanvar_model$omega_bar <- fit$meanvar_model$omega_bar / scale^2
      fit$meanvar_model$mu_lambda_bar <- fit$meanvar_model$mu_lambda_bar * scale
      fit$meanvar_model$mu2_lambda_bar <- fit$meanvar_model$mu2_lambda_bar * scale^2
    }
  }

  return(fit)
}
