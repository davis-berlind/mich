# Expected Squared-Mean Signal

Given that \\E\[y_t\|b,\tau\]=\mu_t=bI(t\geq \tau)\\,
\\P(\tau=t)=\pi_t\\, \\\bar{b}\_t=E\[b\|\tau=t\]\\, and
\\\bar{\omega}\_t=V(b\|\tau=t)\\, `mu2_bar_fn()` calculates
\\E\[\mu^2_t\]\\ as \\\Sigma\_{t'=1}^t (\bar{b}^2\_{t'} + 1 /
\bar{\omega}\_{t'})\pi\_{t'}\\.

## Usage

``` r
mu2_bar_fn(b, omega, prob)
```

## Arguments

- b:

  A numeric vector. Length \\T\\ vector of conditional mean parameters.

- omega:

  A numeric vector. Length \\T\\ vector of conditional variance
  parameters.

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric vector. A length \\T\\ vector of \\E\[\mu_t\]\\.
