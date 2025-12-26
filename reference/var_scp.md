# Variance Single Change-Point Model

Implementation of the Var-SCP model from Berlind, Cappello, and Madrid
Padilla (2025). The function `var_scp()` takes a length \\T\\ vector
\\y\_{1:T}\\ with a single variance change and returns the posterior
distribution of the change-point.

## Usage

``` r
var_scp(y, omega, u_bar, lgamma_u_bar, v, log_pi)
```

## Arguments

- y:

  A numeric vector. \\T\\ observations with a single variance change.

- omega:

  A numeric vector. Known trend component of precision of `y`.

- u_bar:

  A numeric vector. Posterior shape parameters equal to \\u_0 + T - t +
  1\\ for each \\t\\.

- lgamma_u_bar:

  A numeric vector. Log gamma function evaluated at u_bar.

- v:

  A numeric vector. Vector of prior rate parameters for each \\t\\.

- log_pi:

  A numeric vector. Vector of log prior probabilities for the location
  of the change-point.

## Value

A list. A list of posterior parameters including the rate `v_bar`, and
posterior probabilities of the change-point location `pi_bar`.
