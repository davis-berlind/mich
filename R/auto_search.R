multi_auto_search <- function(
    y,
    fit,
    log_pi_l,
    omega,
    trace_omega_omega_bar_inv,
    log_det_omega_omega_bar,
    mean_weights,
    sandwich_weights,
    L_max,
    L_auto,
    fit_intercept,
    fit_scale,
    max_iter,
    tol,
    verbose,
    merge_prob,
    restart,
    n_search,
    increment
) {
  T <- nrow(y)
  d <- ncol(y)

  L <- fit$L
  post_params <- fit$post_params
  mu_0 <- fit$mu_0
  lambda <- fit$lambda
  log_lambda <- log(lambda)
  Q <- fit$Q

  # if components were merged out use auto procedure to increase to desired L
  merge_flag <- (L < L_max & !L_auto)
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components

  # 0. Initilaize fit_best
  # 1. Increase L, if ELBO increases, set fit to current fit, if ELBO decreases,
  #    decrease counter and go back to 1.
  # 2. If counter == 0:
  # 3. If restart == TRUE, reset components to null model, set restart = FALSE
  #    and go back to 1.
  # 4. If restart == FALSE and fit_best is better than fit, return fit_best,
  #    otherwise merge components if no merges return model, otherwise go back
  #    to 1.

  last_restart <- ifelse(restart, 2, Inf)
  elbo_best <- -Inf
  if (L_auto | merge_flag) {
    counter <- ifelse(merge_flag, L_max-L, n_search) # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # store current value of ELBO
    elbo_new <- elbo

    if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))

    # continue search until n_search exhausted or max components exceeded
    while (L < L_max) {
      if (increment > L_max - L) increment <- L_max - L
      # increment dimension of parameters
      for (i in 1:increment) {
        post_params[[L+i]] <- param_init(1, T, d)[[1]]
      }
      if (L > 0) {
        post_params <- post_params[c((L+1):(L+increment), 1:L)]
        if (L_auto) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T+1, ncol = increment), log_pi_l)
        } else {
          log_pi_l <- log_pi_l[, c((L_max - L + 1):L_max, 1:(L_max - L))]
        }
      }
      L <- L + increment

      merged <- FALSE # flag indicating components have been merged
      merge_dist <- ceiling(log(T) / 2)
      while (!merged) {
        fit_new <- mich_matrix_cpp(
          y, mu_0,
          lambda, log_lambda, Q,
          trace_omega_omega_bar_inv, log_det_omega_omega_bar,
          mean_weights, sandwich_weights,
          fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
          omega, log_pi_l,
          post_params
        )

        elbo_new <- fit_new$elbo[length(fit_new$elbo)]
        if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))
        if (merge_dist == 0 | L <= 1) break

        fit_new <- multi_component_merge(
          fit_new, lambda, mean_weights, merge_prob, merge_dist, drop = FALSE
        )

        merge_dist <- merge_dist - 1

        L <- fit_new$L
        post_params <- fit_new$post_params
        prob_order <- order(
          sapply(1:L, function(l) post_params[[l]][["pi_bar"]][T+1]),
          decreasing = TRUE
        )
        post_params <- post_params[prob_order]
        log_pi_l <- log_pi_l[, prob_order, drop = FALSE]

        merged <- fit_new$merged
      }

      # test if model improved or restart ####
      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- ifelse(merge_flag, L_max-L, n_search)
        if (last_restart < Inf) restart <- TRUE
      } else {
        counter <- counter - 1
      }

      if (counter == 0 | L == L_max) {
        if (restart & L > last_restart & !merge_flag) {
          if (verbose) print(paste0("Restarting at L = ", L))
          restart <- FALSE
          last_restart <- L
          post_params <- param_init(L-increment, T, d)
        } else {
          if (elbo > elbo_best) {
            fit_best <- fit
            elbo_best <- elbo
          } else {
            return(fit_best)
          }

          # set fit to best model after current search
          L <- fit$L
          post_params <- fit$post_params
          mu_0 <- fit$mu_0
          lambda <- fit$lambda
          Q <- fit$Q
          if (L_auto) log_pi_l <- log_pi_l[,1:max(1, L), drop=FALSE]

          merged <- FALSE
          merge_dist <- 2
          while (!merged) {
            if (verbose) print("Merging")
            fit_merge <- mich_matrix_cpp(
              y, mu_0,
              lambda, log_lambda, Q,
              trace_omega_omega_bar_inv, log_det_omega_omega_bar,
              mean_weights, sandwich_weights,
              fit_intercept, fit_scale, max_iter, tol, verbose = verbose & !L_auto,
              omega, log_pi_l,
              post_params
            )

            elbo_merge <- fit_merge$elbo[length(fit_merge$elbo)]
            if (elbo_merge > elbo_best) {
              fit_best <- fit_merge
              elbo_best <- elbo_merge
            }
            if (verbose)  print(paste0("L = ", fit_merge$L,"; ELBO = ", elbo_merge))
            if (merge_dist == 0 | L <= 1) break

            fit_merge <- multi_component_merge(
              fit_merge, lambda, mean_weights, merge_prob, merge_dist, drop = TRUE
            )

            merge_dist <- merge_dist - 1

            L <- fit_merge$L
            post_params <- fit_merge$post_params
            prob_order <- order(
              sapply(1:L, function(l) post_params[[l]][["pi_bar"]][T+1]),
              decreasing = TRUE
            )
            post_params <- post_params[prob_order]
            if (L_auto) log_pi_l <- log_pi_l[,1:max(L,1), drop = FALSE]
            else log_pi_l <- log_pi_l[, prob_order, drop = FALSE]

            merged <- fit_merge$merged
          }
          if (merge_dist == 1) {
            return(fit_best)
          } else {
            elbo <- elbo_merge
            counter <- ifelse(merge_flag, L_max-L, n_search)
          }
        }
      }
    }
  }
  return(fit)
}
