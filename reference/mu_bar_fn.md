# Expected Mean Signal

Given that \\E\[y_t\|b,\tau\]=\mu_t=bI(t\geq \tau)\\,
\\P(\tau=t)=\pi_t\\, and \\\bar{b}\_t=E\[b\|\tau=t\]\\, `mu_bar_fn()`
calculates \\E\[\mu_t\]\\ as \\\Sigma\_{t'=1}^t
\bar{b}\_{t'}\pi\_{t'}\\.

## Usage

``` r
mu_bar_fn(b, prob)
```

## Arguments

- b:

  A numeric vector. Length \\T\\ vector of conditional mean parameters.

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric vector. A length \\T\\ vector of \\E\[\mu_t\]\\.
