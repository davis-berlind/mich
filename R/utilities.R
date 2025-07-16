#' Reverse cumulative sum
#'
#' @description
#' `revcumsum()` takes a length \eqn{n} vector \eqn{\{x_i\}_{i=1}^n} and returns \eqn{\{\Sigma_{j=i}^n x_j\}_{i=1}^n}.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector. Reverse cumulative sum of the elements of `x`.
#' @export
#'
#' @examples
#' revcumsum(1:10)
revcumsum <- function(x){
  T <- length(x)
  s <- c(0, cumsum(x[-T]))
  return(x[T] + s[T] - s)
}

#' Probability vector check.
#'
#' @description
#' `prob_check()` throws an error if `probs` is not either one of the string
#' arguments 'uniform' or 'weighted', a length \eqn{T} vector with elements that
#' sum to one, or an \eqn{T \times n} matrix with columns that sum to one.
#'
#' @param probs A character, numeric vector, or numeric matrix. If character,
#'   `probs` should be equal to 'weighted' or 'uniform'. If vector, elements of
#'   `probs` must sum to one. If matrix, rows of `probs` must sum to one.
#' @param n An integer. Number of columns if `probs` is a matrix.
#' @param T An integer. Length or number of rows in `probs`.
#' @return A character, numeric vector, or numeric matrix. If no error is thrown
#'   during evaluation, then `prob_check()` returns the `probs` argument.
#'
prob_check <- function(probs, n, T) {
  if (length(probs) == 1) {
    if (probs %in% c("uniform", "weighted")) return(probs)
  } else if (is.vector(probs) & is.numeric(probs)) {
    if (all(!is.na(probs)) & (length(probs) == T) & is.numeric(probs)) {
      if (round(sum(probs), 10) == 1) return(sapply(1:max(1, n), function(i) probs))
    }
  } else if (is.array(probs) & is.numeric(probs)) {
    if (all(!is.na(probs)) & (nrow(probs) == T & ncol(probs) == n)) {
      if (all(round(colSums(probs), 10) == 1))  return(probs)
    }
  }
  stop(paste0(deparse(substitute(probs)),
              " must be 'uniform', 'weighted', or a length T vector or a T x ",
              deparse(substitute(n)),
              " matrix with columns that sum to one."))
}

#' Logical value check.
#'
#' @description
#' `logical_check()` throws an error if `x` is not equal to either `TRUE` or
#' `FALSE`.
#'
#' @param x A logical.
#' @return A logical. If no error is thrown during evaluation, then
#'   `logical_check()` returns the `x` argument.
#'
logical_check <- function(x) {
  if (is.logical(x) & length(x) == 1) {
    if (!is.na(x)) return(x)
  }
  stop(paste0(deparse(substitute(x)), " must be TRUE or FALSE."))
}

#' Scalar value check.
#'
#' @description
#' `scalar_check()` throws an error if `x` is not equal to a scalar value.
#'
#' @param x A scalar.
#' @return A scalar If no error is thrown during evaluation, then
#'   `scalar_check()` returns the `x` argument.
#'
scalar_check <- function(x) {
  if (is.numeric(x) & length(x) == 1) {
    if (!is.na(x)) {
      if (x > 0) return(x)
    }
  }
  stop(paste0(deparse(substitute(x)), " must be an scalar > 0."))
}

#' Integer value check.
#'
#' @description
#' `integer_check()` throws an error if `n` is not equal to a positive valued
#' integer.
#'
#' @param n An integer.
#' @return An integer. If no error is thrown during evaluation, then
#'   `integer_check()` returns the `n` argument.
#'
integer_check <- function(n) {
  if (is.numeric(n) & length(n) == 1) {
    if (!is.na(n)) {
      if (n == round(n) & n >= 0) return(n)
    }
  }
  stop(paste0(deparse(substitute(n)), " must be an integer >= 0."))
}
