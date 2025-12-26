# MICH Posterior Credible Sets

The function `mich_sets()` takes a T x N matrix of posterior
change-point location probabilities `probs`, a coverage level `level`,
and a max set length, and for column `i` of `probs` returns the MAP
estimator `which.max(probs[,i])` and the smallest set of indices `s_i`
such that `probs[s_i,i] > level` if `length(s_i) < max_length`.

## Usage

``` r
mich_sets(probs, max_length = log(nrow(probs))^1.5, level = 0.9)
```

## Arguments

- probs:

  A numeric Matrix. A T x N matrix of posterior probabilities for the
  location of the change-points.

- max_length:

  A positive scalar. Detection threshold, if a credible set contains
  more that `max_length` indices, then no change is detected. Set equal
  to `log(T)^1.5` by default (see Section 2.5 of Berlind, Cappello, and
  Madrid Padilla (2025)).

- level:

  A scalar. A single number in (0,1) that gives the lower bound for the
  probability that each credible set contains a change-point.

## Value

A list. MAP estimator of each change-point and corresponding credible
set.
