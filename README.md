
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mich <a href="https://davis-berlind.github.io/mich/"><img src="man/figures/logo.png" align="right" height="100" alt="mich website" /></a>

## Fast Bayesian Inference for Change-Point Detection

<!-- badges: start -->
<!-- badges: end -->

`{mich}` is an R package that implements the Multiple Independent
Change-Point (MICH) method introduced in [Berlind, Cappello, Madrid
Padilla (2025)](https://arxiv.org/abs/2507.01558). The package’s main
function `mich()` implements a backfitting procedure to identify changes
in the mean and/or variance of a length $T$ sequence of observations
$\mathbf{y}_{1:T}$. The `mich.fit` object returned by `mich()` provides
a variational approximation to the posterior distribution of the
change-points.

## Installation

The development version of `{mich}` can be installed from GitHub using:

``` r
# install.packages("devtools")
devtools::install_github("davis-berlind/mich")
```
