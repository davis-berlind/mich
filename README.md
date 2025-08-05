# Fast Bayesian Inference for Change-Point Detection

`mich` is an R package that implements the Multiple Independent Change-Point (MICH) method introduced in 
[Berlind, Cappello, Madrid Padilla (2025)](https://arxiv.org/abs/2507.01558). The main function in the package 
`mich()` takes a length $T$ sequence of observations $\mathbf{y}_{1:T}$ with potentially many 
change-points in the mean and variance, and deploys a backfitting procedure to find a variational 
approximation to the posterior distribution of the change-points.

## Installation

You can install the development version of `mich` from GitHub using:

```{r}
# install.packages("devtools")
devtools::install_github("davis-berlind/mich")
```
