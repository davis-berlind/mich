multi_auto_search <- function(
    y,
    fit,
    omega_l,
    log_pi_l,
    L_max,
    L_auto,
    fit_intercept,
    fit_scale,
    max_iter,
    tol,
    verbose,
    merge_level,
    merge_prob,
    restart,
    n_search,
    increment
) {
  T <- nrow(y)
  d <- ncol(y)

  # max times to merge
  merge_counter <- log(T) %/% 2

  lambda <- fit$lambda
  Q <- fit$Q
  mu_0 <- fit$mu_0
  post_params <- fit$post_params
  L <- fit$L

  # log params
  log_lambda <- log(lambda)
  log_omega_l <- log(omega_l)
  d_log_omega_l <- d * log_omega_l

  # eigen value parameters
  eigen_vals <- omega_l + outer((T-1:(T+1)+1), lambda)
  log_prob_weights <- log_pi_l - 0.5 * rowSums(log(eigen_vals))
  Lambda_bar_log_det <- rowSums(log(eigen_vals))
  inv_Lambda_bar_trace <- omega_l / rowSums(eigen_vals)
  mean_weights <- matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE) / eigen_vals
  sandwich_weights <- mean_weights * matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE)

  # if components were merged out use auto procedure to increase to desired L
  merge_flag <- (L < L_max & !L_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components

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
    counter <- ifelse(merge_flag, L_max-L, n_search) # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # store current value of ELBO
    elbo_new <- elbo

    if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))

    # continue search until n_search exhausted or max components exceeded
    while (L < L_max) {
      # increment dimension of parameters
      for (i in 1:increment) {
        post_params[[L+i]] <- param_init(1, T, d)[[1]]
      }
      if (L > 0) {
        post_params <- post_params[c((L+1):(L+increment), 1:L)]
        if (L_auto) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T+1, ncol = increment), log_pi_l)
          log_prob_weights <- cbind(matrix(log_prob_weights[,1], nrow = T+1, ncol = increment), log_prob_weights)
        } else {
          log_pi_l <- log_pi_l[, c((L_max - L + 1):L_max, 1:(L_max - L))]
          log_prob_weights <- log_prob_weights[, c((L_max - L + 1):L_max, 1:(L_max - L))]
        }
      }
      L <- L + increment

      # fit incremented model
      fit_new <- mich_matrix_cpp(
        y, mu_0,
        lambda, log_lambda, Q,
        Lambda_bar_log_det, inv_Lambda_bar_trace,
        mean_weights, sandwich_weights, log_prob_weights,
        fit_intercept, fit_scale, max_iter, tol, verbose = FALSE,
        omega_l, log_omega_l, d_log_omega_l, log_pi_l,
        post_params
      )

      # test if model improved or merge/restart ####
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new,"; Counter: ", counter))

      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- ifelse(merge_flag, L_max-L, n_search)
        if (last_restart < Inf) restart <- TRUE
      } else {
        counter <- counter - 1
      }

      if (counter == 0 | L == L_max) {
        if (restart & !merge_flag) {
          if (verbose) print(paste0("Restarting at L = ", L))
          restart <- FALSE
          last_restart <- L
          post_params <- param_init(L, T, d)
          counter <- ifelse(merge_flag, L_max-L, n_search)
        } else {
          # merging components
          merge_counter <- merge_counter - 1
          no_merges <- TRUE # flag indicating no components were merged
          merged <- FALSE  # flag indicating components have been merged
          if (verbose) print(paste0("Merging. Merge Counter: ", merge_counter))

          # set fit to best model so far
          L <- fit$L
          post_params <- fit$post_params
          mu_0 <- fit$mu_0
          lambda <- fit$lambda
          Q <- fit$Q
          if (L_auto) log_pi_l <- log_pi_l[,1:max(1, L), drop=FALSE]

          eigen_vals <- omega_l + outer((T-1:(T+1)+1), lambda)
          log_prob_weights <- log_pi_l - 0.5 * rowSums(log(eigen_vals))
          Lambda_bar_log_det <- rowSums(log(eigen_vals))
          inv_Lambda_bar_trace <- omega_l / rowSums(eigen_vals)
          mean_weights <- matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE) / eigen_vals
          sandwich_weights <- mean_weights * matrix(lambda, nrow = T+1, ncol = d, byrow = TRUE)

          while (!merged) {
            if (merge_flag & merge_counter == 0) break
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

            if (fit$L < L) {
              no_merges <- FALSE
              L <- fit$L
            }
            post_params <- fit$post_params

            if (L_auto) {
              log_pi_l <- log_pi_l[,1:max(1,L), drop=FALSE]
              log_prob_weights <- log_prob_weights[,1:max(1,L), drop=FALSE]
            }

            merged <- fit$merged
            if (verbose & !merged) print(paste0("Merging to L = ", L))
          }

          elbo <- max(fit$elbo)
          if (elbo > merge_elbo | (merge_flag & merge_counter == 0)) {
            merge_elbo <- elbo
            fit_merge <- fit
          }

          counter <- ifelse(merge_flag, L_max-L, n_search) # reset counter in case searching continues

          # return model if no components to merge or reached max number of merge attempts
          if (no_merges | merge_counter == 0) {
            fit <- fit_merge
            break
          }
        }
      }
    }
  }
  return(fit)
}
