# Multivariate Expected Mean Signal

Given that \\E\[\mathbf{y}\_t\|\mathbf{b},\tau\] = \mathbf{b}I(t\geq
\tau)\\, \\P(\tau = t)=\pi_t\\, and
\\\bar{\mathbf{b}}\_t=E\[\mathbf{b}\|\tau=t\]\\, `multi_mu_bar_fn()`
calculates \\E\[\mathbf{b}I(t\geq \tau)\]\\ as \\\Sigma\_{t'=1}^t
\bar{\mathbf{b}}\_{t'}\pi\_{t'}\\.

## Usage

``` r
multi_mu_bar_fn(b, prob)
```

## Arguments

- b:

  A numeric matrix. \\T\times d\\ matrix of conditional mean parameters.

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric matrix. A \\T\times d\\ matrix of \\E\[\mathbf{b}I(t\geq
\tau)\]\\.
