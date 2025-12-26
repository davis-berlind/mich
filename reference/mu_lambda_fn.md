# Expected Mean-Precision Signal

Given that \\E\[y_t\|b,\tau\]=\mu_t=bI(t\geq \tau)\\,
\\\bar{b}\_t=E\[b\|\tau=t\]\\, \\\bar{\omega}\_t=V(b\|\tau=t)\\,
\\V(y_t\|s,\tau)=1/\lambda_t=1/s^{I(t\geq \tau)}\\,
\\\bar{u}\_t/\bar{v}\_t=E\[s\|\tau=t\]\\, and \\P(\tau=t)=\pi_t\\,
`mu_lambda_fn()` calculates \\E\[\mu_t\lambda_t\]\\ as
\\\Sigma\_{t'=1}^t \bar{b}\_{t'}\bar{v}\_{t'} \pi\_{t'}/
\bar{u}\_{t'}\\.

## Usage

``` r
mu_lambda_fn(b, u, v, prob)
```

## Arguments

- b:

  A numeric vector. Length \\T\\ vector of conditional mean parameters.

- u:

  A numeric vector. Length \\T\\ vector of conditional shape parameters

- v:

  A numeric vector. Length \\T\\ vector of conditional rate parameters

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric vector. A length \\T\\ vector of \\E\[\mu_t\lambda_t\]\\.
