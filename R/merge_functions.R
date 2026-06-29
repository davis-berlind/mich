component_merge <- function(
    fit, L_auto, K_auto, J_auto,
    merge_level, merge_prob,
    omega_l, log_pi_l,
    v_k, u_bar_k, lgamma_u_bar_k, log_pi_k,
    omega_j, v_j, u_bar_j, lgamma_u_bar_j, log_pi_j
){
  T <- length(fit$y)
  merged <- TRUE
  L <- fit$L
  J <- fit$J
  K <- fit$K
  # identify components to merge
  if (L > 1) {
    pi_bar_l <- fit$mean_model$pi_bar
    # only merge columns with credible sets with length less than detect
    keep <- pi_bar_l[T+1,] <= 1 - merge_level
    keep_mat <- matrix(FALSE, ncol = L, nrow = L)
    keep_mat[keep, keep] <- TRUE
    diag(keep_mat) <- FALSE

    # compute pairwise merge probabilities
    merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
    diag(merge_prob_mat) <- 0

    b_bar_l <- fit$mean_model$b_bar
    log_pi_bar_l <- fit$mean_model$log_pi_bar
    omega_bar_l <- fit$mean_model$omega_bar
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

    fit$L <- L
    fit$mean_model$b_bar <- b_bar_l
    fit$mean_model$pi_bar <- pi_bar_l
    fit$mean_model$log_pi_bar <- log_pi_bar_l
    fit$mean_model$omega_bar <- omega_bar_l
    fit$log_pi_l <- log_pi_l
  }

  if (K > 1) {
    pi_bar_k <- fit$var_model$pi_bar

    # only merge columns with credible sets with length less than detect
    keep <- pi_bar_k[T+1,] <= 1 - merge_level
    keep_mat <- matrix(FALSE, ncol = K, nrow = K)
    keep_mat[keep, keep] <- TRUE
    diag(keep_mat) <- FALSE

    # compute pairwise merge probabilities
    merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
    diag(merge_prob_mat) <- 0

    lambda_bar_k <- fit$var_model$lambda_bar
    log_pi_bar_k <- fit$var_model$log_pi_bar
    v_bar_k <- fit$var_model$v_bar

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
      log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda * merge_delta))

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

    fit$K <- K
    fit$var_model$v_bar_k <- v_bar_k
    fit$var_model$pi_bar <- pi_bar_k
    fit$var_model$log_pi_bar <- log_pi_bar_k
    fit$log_pi_k <- log_pi_k
  }

  if (J > 1) {
    pi_bar_j <- fit$meanvar_model$pi_bar

    # only merge columns with credible sets with length less than detect
    keep <- pi_bar_j[T+1,] <= 1 - merge_level
    keep_mat <- matrix(FALSE, ncol = J, nrow = J)
    keep_mat[keep, keep] <- TRUE
    diag(keep_mat) <- FALSE

    # compute pairwise merge probabilities
    merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
    diag(merge_prob_mat) <- 0

    mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
    mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
    lambda_bar_j <- fit$meanvar_model$lambda_bar
    b_bar_j <- fit$meanvar_model$b_bar
    log_pi_bar_j <- fit$meanvar_model$log_pi_bar
    omega_bar_j <- fit$meanvar_model$omega_bar
    v_bar_j <- fit$meanvar_model$v_bar

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
      log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda * merge_delta))

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

    fit$J <- J
    fit$meanvar_model$b_bar <- b_bar_j
    fit$meanvar_model$omega_bar <- omega_bar_j
    fit$meanvar_model$v_bar <- v_bar_j
    fit$meanvar_model$pi_bar <- pi_bar_j
    fit$meanvar_model$log_pi_bar <- log_pi_bar_j
    fit$log_pi_j <- log_pi_j
  }
  fit$merged <- merged
  return(fit)
}

multi_component_merge <- function(
  fit, lambda, mean_weights, merge_level, merge_prob
) {
  T <- nrow(fit$QTr)
  post_params <- fit$post_params
  L <- length(post_params)
  merged <- TRUE
  # identify components to merge
  if (L > 1) {
    # extract posterior change-point probabilities
    pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

    # only merge columns with credible sets with pi_bar_T+1 < merge_level
    keep <- pi_bar_l[T+1,] <= 1 - merge_level
    keep_mat <- matrix(FALSE, ncol = L, nrow = L)
    keep_mat[keep, keep] <- TRUE
    diag(keep_mat) <- FALSE

    # compute pairwise merge probabilities
    merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
    diag(merge_prob_mat) <- 0

    # merge components with pairwise merge probabilities > merge_prob
    while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
      merged <- FALSE
      L <- L - 1

      # identify components with largest pairwise merge probabilities
      merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
      if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
        merge_dex <- merge_dex[2:1]
      }

      # store merged mean parameters
      post_params[[merge_dex[1]]][["QTb_bar"]] <- post_params[[merge_dex[1]]][["QTb_bar"]] + post_params[[merge_dex[2]]][["QTb_bar"]]
      post_params[[merge_dex[1]]][["b_bar"]] <- post_params[[merge_dex[1]]][["b_bar"]] + post_params[[merge_dex[2]]][["b_bar"]]
      post_params[[merge_dex[1]]][["QTmu_bar"]] <- multi_mu_bar_fn(
        post_params[[merge_dex[1]]][["QTb_bar"]],
        post_params[[merge_dex[1]]][["pi_bar"]]
      )
      post_params[[merge_dex[1]]][["muTLmu_bar_l"]] <- muTLmu_bar_fn(
        post_params[[merge_dex[1]]][["QTb_bar"]],
        post_params[[merge_dex[1]]][["pi_bar"]],
        lambda,
        mean_weights
      )

      # subtract out mean signals of merged components
      # drop merged component
      post_params[merge_dex[2]] <- NULL
      pi_bar_l <- pi_bar_l[,-merge_dex[2]]
      merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
      keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
    }
  }

  fit$L <- L
  fit$post_params <- post_params
  fit$merged <- merged
  return(fit)
}
