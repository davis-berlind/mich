# Single Change-Point Posterior Credible Set

The function `cred_set()` takes a length T vector of posterior
change-point location probabilities `prob` and a coverage level `level`,
and returns the smallest set of indices `s` such that `prob[s] > level`.

## Usage

``` r
cred_set(prob, level)
```

## Arguments

- prob:

  A numeric vector. A vector of posterior probabilities for the location
  of the change-point.

- level:

  A scalar. A single number in (0,1) that gives the lower bound for the
  probability that the credible set contains a change-point.

## Value

A vector. A Level `level` posterior credible set for the location of a
single change-point.
