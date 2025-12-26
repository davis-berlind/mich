# Log Weighted Mean-SCP Prior

Log weighted prior for the Mean-SCP model as described in Appendix C.2
of Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that
in the absence of a change-point, the posterior probabilities are
approximately uniform.

## Usage

``` r
log_mean_prior(T, d)
```

## Arguments

- T:

  An integer. Number of observations.

- d:

  An integer. Dimension each observation.

## Value

A numeric vector. Log prior probabilities.
