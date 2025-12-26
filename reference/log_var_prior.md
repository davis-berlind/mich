# Log Weighted Var-SCP Prior

Log weighted prior for the Var-SCP model as described in Appendix C.2 of
Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that in
the absence of a change-point, the posterior probabilities are
approximately uniform.

## Usage

``` r
log_var_prior(T)
```

## Arguments

- T:

  An integer. Number of observations.

## Value

A numeric vector. Log prior probabilities.
