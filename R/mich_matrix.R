#' Multivariate Multiple Independent Change-Point (MICH) Model
#'
#' @description
#' Fits the multivariate version of the MICH model for mean change-points.
#' Number of change-points can either be fixed or `mich_matrix()` will search
#' for the number of changes that maximizes the ELBO when `L_auto == TRUE`.
#'
#' @param y A numeric matrix. \eqn{T \times d} matrix of observations.
#' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
#'   intercept is estimated, otherwise it is assumed  that
#'   \eqn{\boldsymbol{\mu}_0 = \mathbf{0}}.
#' @param fit_scale A logical. If `fit_scale == TRUE`, then the precision matrix
#'   \eqn{\Lambda = E[\mathbf{y}_t\mathbf{y}'_t]^{-1}} is estimated using
#'   \eqn{\hat{\Lambda}^{-1} = \frac{1}{2(T-1)} \Sigma_{t=1}^{T-1} (\mathbf{y}_{t+1} - \mathbf{y}_{t})(\mathbf{y}_{t+1} - \mathbf{y}_{t})'},
#'   otherwise it is assumed that \eqn{\Lambda = \mathbf{I}_d}.
#' @param standardize A logical. If `standardize == TRUE`, then `y` is centered
#'    and rescaled before fitting.
#' @param L An integer. Number of mean-scp components included in model. If
#'   `L_auto == TRUE` then `L` lower bounds the number of change-points in the
#'   model.
#' @param L_auto A logical. If `L_auto == TRUE`, then `mich_matrix()` returns
#'   the \eqn{L} between `L` and `L_max` that maximizes the ELBO (see Appendix
#'   C.4 of Berlind, Cappello, Madrid Padilla (2025)).
#' @param L_max L An integer. If `L_auto == TRUE` then `L_max` upper bounds the
#'   number of change-points included in the model.
#' @param pi_l_weighted A logical. If `pi_l_weighted == TRUE`, then the weighted
#'  priors specified in Appendix C.2 of Berlind, Cappello, Madrid Padilla (2025)
#'  are used.
#' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
#' @param verbose A logical. If `verbose == TRUE` and `L_auto == FALSE`, then
#'   the value of the ELBO is printed every 5000th iteration. If
#'   `verbose == TRUE` and `L_auto == TRUE`, the the value of the ELBO is
#'   printed for each \eqn{L} as `mich_matrix()` searches over \[`L`, `L_max`\].
#' @param max_iter An integer. Maximum number of iterations. If ELBO does not
#'   converge before `max_iter` is reached, then `converged == FALSE` in the
#'   returned fit object.
#' @param reverse A logical. If `reverse == TRUE` then MICH is fit to
#'   \eqn{\mathbf{y}_{T:1}} and the model parameters are reversed in
#'    post-processing.
#' @param detect A scalar. The detection criteria. The \eqn{\ell^{\text{th}}}
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
#' @param log_pi_l A numeric matrix. A \eqn{T \times L} matrix of prior log
#'   change-point location probabilities for each of the \eqn{L} mean
#'   change-points.
#'
#' @return A list. Parameters of the variational approximation the MICH
#' posterior distribution, including:
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
mich_matrix <- function(y, fit_intercept, fit_scale, standardize,
                        L, L_auto, L_max, pi_l_weighted,
                        tol, verbose, max_iter, reverse,
                        detect, merge_level, merge_prob,
                        restart, n_search, increment,
                        omega_l, log_pi_l) {

  #### set up ####
  # Flag for new model
  refit <- FALSE

  # store original y
  y_raw <- y

  # calculate dimensions of y
  T <- nrow(y)
  d <- ncol(y)

  # min prob to keep component when restarting
  keep_level <- 0.9

  # max times to merge
  merge_counter = log(T) %/% 2

  #### standardize data ####
  if (fit_scale | standardize) {
    # estimate precision matrix
    y_diff <- diff(y)
    y_diff_norm <- sqrt(rowSums(y_diff^2))

    # remove outliers due to big mean changes
    y_diff <- y_diff[y_diff_norm <= stats::quantile(y_diff_norm, p = 0.75) +  1.5 * stats::IQR(y_diff_norm), ]

    # estimate variance
    Sigma <- stats::var(y_diff) / 2
    Sigma_eigen <- eigen(Sigma)
    if (any(Sigma_eigen$values <= 0)) {
      warning("Var(y) is singular. Consider removing collinear columns.")
      Sigma_eigen$values <- Sigma_eigen$values + 1e-5
    }
    Lambda_sqrt <- Sigma_eigen$vectors %*% diag(1 / sqrt(Sigma_eigen$values)) %*% t(Sigma_eigen$vectors)

    # center data
    if (standardize) {
      center <-  apply(y, 2, stats::median)
      y <- y - center
    }

    # rescale data
    y <- y %*% Lambda_sqrt
  }

  #### initialize mu_0 ####
  if (fit_intercept) {
    mu_0 <- colMeans(y[1:ceiling(log(T)),,drop=FALSE])
  } else mu_0 <- rep(0.0, d)

  #### initializing posterior parameters ####
  omega_bar_l <- omega_l + T - 1:T + 1
  log_omega_bar_l <- log(omega_bar_l)
  post_params <- list()
  if (L > 0) {
    for (l in 1:L) {
      post_params[[l]] <- list(pi_bar = rep(1 / T, T),
                               log_pi_bar = rep(0.0, T),
                               b_bar = matrix(0.0, nrow = T, ncol = d))
    }
  }

  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  while (!merged) {
    fit <- multi_mich_cpp(y, mu_0, fit_intercept, refit,
                          max_iter, tol, verbose = verbose & !L_auto,
                          omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                          post_params)

    merged <- TRUE
    refit <- TRUE

    # identify components to merge
    if (L > 1) {
      # extract posterior change-point probabilities
      pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE

      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0

      # store residual
      merge_residual <- fit$residual

      # merge components with pairwise merge probabilities > merge_prob
      while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        L <- L - 1

        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        # add back mean signals of merged components
        mu_bar_merge <- Reduce("+", lapply(merge_dex, function(l) multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])))
        merge_residual <- merge_residual + mu_bar_merge

        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
          pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
          post_params[[merge_dex[1]]][["pi_bar"]] <- post_params[[merge_dex[2]]][["pi_bar"]]
          post_params[[merge_dex[1]]][["log_pi_bar"]] <- post_params[[merge_dex[2]]][["log_pi_bar"]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }

        # store merged mean parameters
        post_params[[merge_dex[1]]][["b_bar"]] <- post_params[[merge_dex[1]]][["b_bar"]] + post_params[[merge_dex[2]]][["b_bar"]]

        # subtract out mean signals of merged components
        merge_residual <- merge_residual - multi_mu_bar_fn(post_params[[merge_dex[1]]][["b_bar"]], post_params[[merge_dex[1]]][["pi_bar"]])

        # drop merged component
        post_params[merge_dex[2]] <- NULL
        pi_bar_l <- pi_bar_l[,-merge_dex[2]]
        merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        if (L_auto) {
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }
      }
    }
    if (verbose & !merged) print(paste0("Merging to L = ", L))
  }

  # if components were merged out use auto procedure to increase to desired L
  merge_flag <- (L < L_max & !L_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components

  #### auto procedure ####

  # 1. Increase L, if ELBO increases, set restart == TRUE and increase L again
  # 2. If ELBO decreases, decrease counter and fit another component
  # 3. If ELBO decreases, counter == 0, and restart == TRUE, reset components
  #    to null model (except for this with concentrated probs) and refit
  # 4. If ELBO decreases, counter == 0, restart == FALSE, and merge_counter > 0,
  #    set fit to best model so far and merge components if no merges return
  #    model, otherwise decrease merge_count and go back to 1
  # 5. If merge_count = 0 return fit

  last_restart <- ifelse(restart, 2, Inf)

  if (L_auto | merge_flag) {
    refit <- (L > 0)
    counter <- n_search # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # store current value of ELBO
    elbo_new <- elbo

    if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))

    # continue search until n_search exhausted or max components exceeded
    while (L < L_max) {
      # increment dimension of parameters
      for (i in 1:increment) {
        post_params[[L+i]] <- list(pi_bar = rep(1 / T, T),
                                   log_pi_bar = rep(0.0, T),
                                   b_bar = matrix(0.0, nrow = T, ncol = d))
        if (L > 1) post_params <- post_params[c(L+i, 1:(L+i-1))]
      }

      L <- L + increment
      if (L > 1 & L_auto) {
        log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
      }

      # fit incremented model
      fit_new <- multi_mich_cpp(y, mu_0, fit_intercept, refit,
                                max_iter, tol, verbose = FALSE,
                                omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                                post_params)

      # test if model improved or merge/restart ####
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new,"; Counter: ", counter))

      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE
      } else if (counter == 0 & restart) {
        if (verbose) print(paste0("Restarting at L = ", L))
        restart <- FALSE
        last_restart <- L
        counter <- n_search

        pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

        # reorder by longest blocks first
        chp <- apply(pi_bar_l, 2, which.max)
        chp_order <- order(chp)
        chp <- chp[chp_order]

        post_params <- post_params[chp_order]
        diff_order <- order(diff(c(chp, T)))
        post_params <- post_params[diff_order]

        # identify components with max prob > keep_level
        keep <- apply(pi_bar_l, 2, max) > keep_level

        # reset components with max prob < keep_level to null components
        if (sum(keep) < L) {
          keep_inc <- max(sum(!keep) - increment, 0)
          L <- sum(keep) + keep_inc
          L <- L - increment

          log_pi_l <- sapply(1:L, function(i) log_pi_l[,1, drop = FALSE])
          post_params <- post_params[c(which(!keep)[seq_len(keep_inc)], which(keep))]
          for (l in seq_len(keep_inc)) {
            post_params[[l]] <- list(pi_bar = rep(1 / T, T),
                                     log_pi_bar = rep(0.0, T),
                                     b_bar = matrix(0.0, nrow = T, ncol = d))
          }
        }
      } else if (counter == 0 | L == L_max) {
        # merging components
        counter <- n_search # reset counter in case searching continues
        merge_counter <- merge_counter - 1
        no_merges <- TRUE # flag indicating no components were merged
        merged <- FALSE  # flag indicating components have been merged
        if (verbose) print(paste0("Merging. Merge Counter: ", merge_counter))

        # set fit to best model so far
        L <- fit$L
        post_params <- fit$post_params

        while (!merged) {
          fit <- multi_mich_cpp(y, mu_0, fit_intercept, refit = TRUE,
                                max_iter, tol, verbose = FALSE,
                                omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                                post_params)

          merged <- TRUE

          # identify components to merge
          if (L > 1) {
            # extract posterior change-point probs
            pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

            # only merge columns with credible sets with length less than detect
            cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = L, nrow = L)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE

            # compute pairwise merge probabilities
            merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
            diag(merge_prob_mat) <- 0

            # store residual
            merge_residual <- fit$residual

            # merge components with pairwise merge probabilities > merge_prob
            while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
              no_merges <- FALSE
              merged <- FALSE
              L <- L - 1

              # identify components with largest pairwise merge probabilities
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

              # add back mean signals of merged components
              mu_bar_merge <- Reduce("+", lapply(merge_dex, function(l) multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])))
              merge_residual <- merge_residual + mu_bar_merge

              # keep probabilities of component with largest posterior probability
              if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
                pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
                post_params[[merge_dex[1]]][["pi_bar"]] <- post_params[[merge_dex[2]]][["pi_bar"]]
                post_params[[merge_dex[1]]][["log_pi_bar"]] <- post_params[[merge_dex[2]]][["log_pi_bar"]]
                merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
              }

              # store merged mean parameters
              post_params[[merge_dex[1]]][["b_bar"]] <- post_params[[merge_dex[1]]][["b_bar"]] + post_params[[merge_dex[2]]][["b_bar"]]

              # subtract out mean signals of merged components
              merge_residual <- merge_residual - multi_mu_bar_fn(post_params[[merge_dex[1]]][["b_bar"]], post_params[[merge_dex[1]]][["pi_bar"]])

              # drop merged component
              post_params[merge_dex[2]] <- NULL
              pi_bar_l <- pi_bar_l[,-merge_dex[2]]
              merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
              keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
              if (L_auto) {
                log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
              }
            }
          }
          if (verbose & !merged) print(paste0("Merging to L = ", L))
        }

        elbo <- max(fit$elbo)
        if (elbo > merge_elbo) {
          merge_elbo <- elbo
          fit_merge <- fit
        }

        # return model if no components to merge or reached max number of merge attempts
        if (no_merges | merge_counter == 0) {
          fit <- fit_merge
          break
        }
      } else {
        counter <- counter - 1
      }
    }
  }

  #### reverse model if reverse == TRUE ####
  if (reverse){
    if(verbose) print("reversing model")
    y_raw <- y_raw[T:1,]
    L <- fit$L

    # reversed residuals, variance, and intercepts
    r_bar <- fit$residual[T:1,]
    mu_0 <- fit$mu[T,]

    # don't reverse weighted priors
    if (!pi_l_weighted) {
      log_pi_l <- log_pi_l[T:1,1:max(1,L),drop = FALSE]
    } else {
      log_pi_l <- sapply(1:max(1,L), function(i) log_pi_l[,1])
    }

    # reversing mean components
    post_params <- fit$post_params
    if (L > 0) pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

    for (l in seq_len(L)) {
      tau_l <- which.max(pi_bar_l[,l])
      mu_bar_l <- multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])
      r_bar_l <- fit$residual + mu_bar_l - matrix(colMeans(mu_bar_l[tau_l:T,,drop = FALSE]), nrow = T, ncol = d, byrow = TRUE)
      r_bar_l <- r_bar_l[T:1,]

      fit_scp <- multi_mean_scp(r_bar_l, omega_bar_l, log_omega_bar_l, log_pi_l[,l])

      pi_bar_l[,l] <- fit_scp$pi_bar
      post_params[[l]][["pi_bar"]] <- fit_scp$pi_bar
      post_params[[l]][["log_pi_bar"]] <- fit_scp$log_pi_bar
      post_params[[l]][["b_bar"]] <- fit_scp$b_bar
    }

    # refit model and merge components
    merged <- FALSE # flag indicating components have been merged
    while (!merged) {
      fit <- multi_mich_cpp(y[T:1,], mu_0, fit_intercept, refit = TRUE,
                            max_iter, tol, verbose = FALSE,
                            omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                            post_params)

      merged <- TRUE

      # identify components to merge
      if (L > 1) {
        # extract posterior change-point probabilities
        pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

        # only merge columns with credible sets with length less than detect
        cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        keep_mat <- matrix(FALSE, ncol = L, nrow = L)
        keep_mat[keep, keep] <- TRUE
        diag(keep_mat) <- FALSE

        # compute pairwise merge probabilities
        merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
        diag(merge_prob_mat) <- 0

        # store residual
        merge_residual <- fit$residual

        # merge components with pairwise merge probabilities > merge_prob
        while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
          merged <- FALSE
          L <- L - 1

          # identify components with largest pairwise merge probabilities
          merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

          # add back mean signals of merged components
          mu_bar_merge <- Reduce("+", lapply(merge_dex, function(l) multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])))
          merge_residual <- merge_residual + mu_bar_merge

          # keep probabilities of component with largest posterior probability
          if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
            pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
            post_params[[merge_dex[1]]][["pi_bar"]] <- post_params[[merge_dex[2]]][["pi_bar"]]
            post_params[[merge_dex[1]]][["log_pi_bar"]] <- post_params[[merge_dex[2]]][["log_pi_bar"]]
            merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
            merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
          }

          # store merged mean parameters
          post_params[[merge_dex[1]]][["b_bar"]] <- post_params[[merge_dex[1]]][["b_bar"]] + post_params[[merge_dex[2]]][["b_bar"]]

          # subtract out mean signals of merged components
          merge_residual <- merge_residual - multi_mu_bar_fn(post_params[[merge_dex[1]]][["b_bar"]], post_params[[merge_dex[1]]][["pi_bar"]])

          # drop merged component
          post_params[merge_dex[2]] <- NULL
          pi_bar_l <- pi_bar_l[,-merge_dex[2]]
          merge_prob_mat<- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
          keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
          if (L_auto) {
            log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
          }
        }
      }
      if (verbose & !merged) print(paste0("Merging to L = ", L))
    }
  }

  #### return model ####
  class(fit) <- "mich.fit"
  fit$y <- y_raw

  # rescale data to original units
  if (fit_scale | standardize) {
    Sigma_sqrt <- Sigma_eigen$vectors %*% diag(sqrt(Sigma_eigen$values)) %*% t(Sigma_eigen$vectors)
    # rescale (and center) posterior parameters
    fit$mu_0 <- c(fit$mu_0 %*% Sigma_sqrt)
    fit$mu <- fit$mu %*% Sigma_sqrt
    if (standardize) {
      fit$mu_0 <- fit$mu_0 + center # recenter mu_0
      fit$mu <- fit$mu + matrix(center, nrow = T, ncol = d, byrow = TRUE)
    }
    for (l in seq_len(fit$L)) {
      fit$post_params[[l]][["b_bar"]] <- post_params[[l]][["b_bar"]]  %*% Sigma_sqrt
    }
    fit$Sigma <- Sigma
  }

  if (fit$L > 0) fit$pi_bar <- sapply(1:fit$L, function(l) fit$post_params[[l]][["pi_bar"]])

  return(fit)
}
