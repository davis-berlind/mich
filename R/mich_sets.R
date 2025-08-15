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
  order(prob, decreasing = TRUE)[1:which.max(cumsum(prob[order(prob, decreasing = TRUE)]) > level)]
}

#' MICH Posterior Credible Sets
#'
#' @description
#' The function `mich_sets()` takes a T x N matrix of posterior change-point
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
mich_sets <- function(probs, max_length = log(nrow(probs))^1.5, level = 0.9) {

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
