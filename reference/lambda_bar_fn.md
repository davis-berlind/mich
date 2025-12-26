# Expected Precision Signal

Given that \\V(y_t\|s,\tau)=1/\lambda_t=1/s^{I(t\geq \tau)}\\,
\\P(\tau=t)=\pi_t\\, and \\\bar{u}\_t/\bar{v}\_t=E\[s\|\tau=t\]\\,
`lambda_bar_fn()` calculates \\E\[\lambda_t\]\\ as \\1 -
\Sigma\_{t'=1}^t \pi\_{t'}(1 - \bar{v}\_{t'} / \bar{u}\_{t'})\\.

## Usage

``` r
lambda_bar_fn(u, v, prob)
```

## Arguments

- u:

  A numeric vector. Length \\T\\ vector of conditional shape parameters

- v:

  A numeric vector. Length \\T\\ vector of conditional rate parameters

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric vector. A length \\T\\ vector of \\E\[\lambda_t\]\\.
