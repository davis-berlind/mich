# Expected Squared-Mean-Precision Signal

Given that \\E\[y_t\|b,\tau\]=\mu_t=bI(t\geq \tau)\\,
\\\bar{b}\_t=E\[b\|\tau=t\]\\, \\V(y_t\|s,\tau)=1/\lambda_t=1/s^{I(t\geq
\tau)}\\, \\\bar{u}\_t/\bar{v}\_t=E\[s\|\tau=t\]\\, and
\\P(\tau=t)=\pi_t\\,
[`mu_lambda_fn()`](https://davis-berlind.github.io/mich/reference/mu_lambda_fn.md)
calculates \\E\[\mu^2_t\lambda_t\]\\ as \\\Sigma\_{t'=1}^t
\pi\_{t'}(\bar{b}^2\_{t'}\bar{v}\_{t'} / \bar{u}\_{t'} +
1/\bar{\omega}\_{t'})\\.

## Usage

``` r
mu2_lambda_fn(b, omega, u, v, prob)
```

## Arguments

- b:

  A numeric vector. Length \\T\\ vector of conditional mean parameters.

- omega:

  A numeric vector. Length \\T\\ vector of conditional variance
  parameters.

- u:

  A numeric vector. Length \\T\\ vector of conditional shape parameters

- v:

  A numeric vector. Length \\T\\ vector of conditional rate parameters

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric vector. A length \\T\\ vector of \\E\[\mu^2_t\lambda_t\]\\.
