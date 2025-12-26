# Probability vector check.

`prob_check()` throws an error if `probs` is not either one of the
string arguments 'uniform' or 'weighted', a length \\T\\ vector with
elements that sum to one, or an \\T \times n\\ matrix with columns that
sum to one.

## Usage

``` r
prob_check(probs, n, T)
```

## Arguments

- probs:

  A character, numeric vector, or numeric matrix. If character, `probs`
  should be equal to 'weighted' or 'uniform'. If vector, elements of
  `probs` must sum to one. If matrix, rows of `probs` must sum to one.

- n:

  An integer. Number of columns if `probs` is a matrix.

- T:

  An integer. Length or number of rows in `probs`.

## Value

A character, numeric vector, or numeric matrix. If no error is thrown
during evaluation, then `prob_check()` returns the `probs` argument.
