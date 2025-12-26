# Mean Single Change-Point Model

Implementation of the univariate Mean-SCP model from Berlind, Cappello,
and Madrid Padilla (2025). The function `mean_scp()` takes a length
\\T\\ vector \\y\_{1:T}\\ with a single mean change and returns the
posterior distribution of the change-point.

## Usage

``` r
mean_scp(y, lambda, omega, log_pi)
```

## Arguments

- y:

  A numeric vector. \\T\\ observations with a single variance change.

- lambda:

  A numeric vector. Precision vector of `y`.

- omega:

  A positive scalar. Prior precision parameter.

- log_pi:

  A numeric vector. Vector of log prior probabilities for the location
  of the change-point.

## Value

A list. A list of posterior parameters including the mean `b_bar`,
precision `omega_bar`, and posterior probabilities of the change-point
location `pi_bar`.
