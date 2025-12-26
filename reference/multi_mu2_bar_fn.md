# Multivariate Expected Squared-Mean Signal

Given that \\E\[\mathbf{y}\_t\|\mathbf{b},\tau\] = \mathbf{b}I(t\geq
\tau)\\, \\P(\tau = t)=\pi_t\\, and for some \\1\leq i \leq d\\,
\\\bar{b}\_{it}=E\[b_i\|\tau=t\]\\ and
\\\bar{\omega}\_{it}=V(b_i\|\tau=t)\\, `multi_mu2_bar_fn()` calculates
\\E\[b_i^2I(t\geq \tau)\]\\ as \\\Sigma\_{t'=1}^t (\bar{b}^2\_{it'} +
1/\bar{\omega}\_{it'}) \pi\_{t'}\\.

## Usage

``` r
multi_mu2_bar_fn(b, omega, prob)
```

## Arguments

- b:

  A numeric matrix. \\T\times d\\ matrix of conditional mean parameters.

- omega:

  A numeric vector. Length \\T\\ vector of conditional variance
  parameters.

- prob:

  A numeric vector. Vector of change-point location probabilities.

## Value

A numeric matrix. A \\T\times d\\ matrix of \\E\[b^2_iI(t\geq \tau)\]\\.
