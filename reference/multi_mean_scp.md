# Mulitvariate Mean Single Change-Point Model

Implementation of the multivariate Mean-SCP model from Berlind,
Cappello, and Madrid Padilla (2025). The function
[`mean_scp()`](https://davis-berlind.github.io/mich/reference/mean_scp.md)
takes a \\T\times d\\ matrix \\\mathbf{y}\_{1:T}\\ with a single mean
change and returns the posterior distribution of the change-point.

## Usage

``` r
multi_mean_scp(y, omega_bar, log_omega_bar, log_pi)
```

## Arguments

- y:

  A numeric matrix. \\T\times d\\ matrix of observations with a single
  mean change.

- omega_bar:

  A numeric vector. Posterior precision parameters \\\bar{\omega}\_t\\
  such that \\\bar{\Omega}\_t = \bar{\omega}\_t\mathbf{I}\_d\\.

- log_omega_bar:

  A numeric vector. Log of `omega_bar`.

- log_pi:

  A numeric vector. Vector of log prior probabilities for the location
  of the change-point.

## Value

A list. A list of posterior parameters including the mean `b_bar`, and
posterior probabilities of the change-point location `pi_bar`.
