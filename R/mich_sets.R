#' Single Change-Point Posterior Credible Set
#'
#' @description
#' The function `cred_set()` takes a length T vector of posterior change-point
#' location probabilities `prob` and a coverage level `level`, and returns the
#' smallest set of indices `s` such that `prob[s] > level`.
#'
#' @param prob  A numeric vector. A vector of posterior probabilities for the
#'   location of the change-point.
#' @param level A scalar. A single number in (0,1) that gives the lower bound
#'   for the probability that the credible set contains a change-point.
#'
#' @return A vector. A Level `level` posterior credible set for the location of
#' a single change-point.
#'
cred_set <- function(prob, level) {
  prob_order <- order(prob, decreasing = TRUE)
  prob_order[1:which.max(cumsum(prob[prob_order]) > level)]
}

#' Marginal Posterior Credible Sets
#'
#' @description
#' The function `marginal_sets()` takes a T x N matrix of posterior change-point
#' location probabilities `probs`, a coverage level `level`, and a max set
#' length, and for column `i` of `probs` returns the MAP estimator
#' `which.max(probs[,i])` and the smallest set of indices `s_i` such that
#' `probs[s_i,i] > level` if `length(s_i) < max_length`.
#'
#' @param probs A numeric Matrix. A T x N matrix of posterior
#'   probabilities for the location of the change-points.
#' @param level A scalar. A single number in (0,1) that gives the lower bound
#'   for the probability that each credible set contains a change-point.
#' @param max_length A positive scalar. Detection threshold, if a credible set
#'   contains more that `max_length` indices, then no change is detected. Set
#'   equal to `log(T)^1.5` by default (see Section 2.5 of Berlind,
#'   Cappello, and Madrid Padilla (2025)).
#'
#' @return A list. MAP estimator of each change-point and corresponding credible
#' set.
#'
#' @export
#'
marginal_sets <- function(
    probs,
    max_length = log(nrow(probs))^1.5,
    level = 0.9
  ) {
  # initialize credible sets and changepoints
  cs <- list()
  est_cp <- numeric(0)

  cred_sets <- apply(probs, 2, cred_set, level = level, simplify = FALSE)

  # drop probs/sets that are longer than max_length
  keep <- which(sapply(cred_sets, length) <= max_length)
  if (length(keep) > 0) {
    est_cp <- apply(probs[, keep, drop = FALSE], 2, which.max)
    cs <- lapply(keep, function(i) cred_sets[[i]])
    # order change-points
    cs <- lapply(order(est_cp), function(i) cs[[i]])
    est_cp <- est_cp[order(est_cp)]
  }
  return(list(cp = est_cp, sets = lapply(cs, function(set) set[order(set)])))
}

#' MICH Posterior Credible Sets
#'
#' @description
#' The function `mich_sets()` takes a T x L matrix of posterior change-point
#' location probabilities `probs`, a coverage level `alpha`, and returns the
#' modified MAP estimator (36) in Berlind et al. (2025). If `detect==TRUE`, the
#' credible sets conditional on detection in (38) are returned, otherwise the
#' unconditional sets (37) are returned.
#'
#' @param probs A numeric Matrix. A T x L matrix of posterior
#'    probabilities for the location of the change-points.
#' @param alpha A scalar. A single number in (0,1) that gives the lower bound
#'    for the probability that each credible set contains a change-point.
#' @param detect A logical. If `TRUE`, the credible sets conditional on
#'    detection are returned.
#'
#' @return A list. MAP estimator of each change-point and corresponding credible
#' set.
#'
#' @export
#'
mich_sets <- function(probs, alpha = 0.95, detect = FALSE) {
  if (is.vector(probs)) probs <- matrix(probs, nrow = length(probs), ncol = 1)
  L <- ncol(probs)
  T <- nrow(probs) - 1

  # determine which columns detect change-points
  detected <- probs[T+1,] <= 0.5

  if (detect) {
    probs <- probs[,detected,drop=FALSE]
    probs <- probs[1:T,,drop=FALSE] / matrix(1 - probs[T+1,], ncol = ncol(probs), nrow = T, byrow = TRUE)
  }

  # calculate sizes of posterior sets
  set_sizes <- set_size(probs, alpha)

  # modified map estimator
  cp <- c()
  sets <- list()
  size_dex <- 0
  for (l in seq_len(L)) {
    if (detected[l] || !detect) {
      size_dex <- size_dex + 1
      order_l <- order(probs[, size_dex], decreasing = TRUE)
      set <- order_l[seq_len(set_sizes[size_dex])]
      if (detected[l]) {
        cp <- c(cp, which.max(probs[1:T, size_dex]))
      } else {
        cp <- c(cp, T + 1)
      }
      sets[[size_dex]] <- set
    }
  }
  if (length(cp) > 1) {
    cp_order <- order(cp)
    cp <- cp[cp_order]
    sets <- sets[cp_order]
  }
  return(list(cp = cp, sets = sets))
}

#' Size of Joint Posterior Credible Sets
#'
#' @description
#' Implementation of the greedy algorithm described in Supplement S3.2 of
#' Berlind et al. (2025).
#'
#' @param probs A numeric Matrix. A T x L matrix of posterior
#'    probabilities for the location of the change-points.
#' @param alpha A scalar. A single number in (0,1) that gives the lower bound
#'    for the probability that each credible set contains a change-point.
#'
#' @return A vector. Sizes of joint alpha-level credible sets.
#'
set_size <- function(probs, alpha) {
  T <- nrow(probs)
  L <- ncol(probs)
  ordered_probs <- matrix(nrow = T, ncol = L)
  ordered_gain <- matrix(nrow = T-1, ncol = L)
  for (l in 1:L) {
    ordered_probs[,l] <- probs[order(probs[,l], decreasing = TRUE),l]
    ordered_gain[,l] <- log(ordered_probs[-1, l]) - log(cumsum(ordered_probs[-T, l]))
  }
  ordered_probs <- rbind(ordered_probs, 0)
  ordered_gain <- rbind(ordered_gain, -Inf)

  cand_prob <- ordered_probs[2,]
  cand_gain <- ordered_gain[1,]
  set_sizes <- rep(1, L)
  comp_total <- ordered_probs[1,]
  total <- prod(comp_total)

  while(total < alpha) {
    l <- which.max(cand_gain)
    set_sizes[l] <- set_sizes[l] + 1
    comp_total[l] <- comp_total[l] + cand_prob[l]
    total <- prod(comp_total)

    cand_prob[l] <- ordered_probs[set_sizes[l]+1, l]
    cand_gain[l] <- ordered_gain[set_sizes[l], l]
  }
  return(set_sizes)
}
