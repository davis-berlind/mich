# Mean-Variance Single Change-Point Model

Implementation of the MeanVar-SCP model from Berlind, Cappello, and
Madrid Padilla (2025). The function `meanvar_scp()` takes a length \\T\\
vector \\y\_{1:T}\\ with a single joint mean and variance change and
returns the posterior distribution of the change-point.

## Usage

``` r
meanvar_scp(y, lambda, omega, u_bar, lgamma_u_bar, v, log_pi)
```

## Arguments

- y:

  A numeric vector. \\T\\ observations with a single joint mean and
  variance change.

- lambda:

  A numeric vector. Known trend component of precision of `y`.

- omega:

  A positive scalar. Prior precision parameter.

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

A list. A list of posterior parameters including the mean `b_bar`,
precision `omega_bar`, rate `v_bar`, and posterior probabilities of the
change-point location `pi_bar`.
