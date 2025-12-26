# Multivariate MICH Algorithm

Implementation of Algorithm 3 from Berlind, Cappello, and Madrid Padilla
(2025). This algorithm takes a sequence of \\d\\-dimensional
observations \\\mathbf{y}\_{1:T}\\, and iteratively fits \\L\\ mean-scp
models resulting in a variational approximation to the posterior
distribution of the \\L\\ mean change-points. The algorithm terminates
once the percentage increase in the ELBO falls below `tol`.

## Usage

``` r
multi_mich_cpp(
  y,
  mu_0,
  fit_intercept,
  refit,
  max_iter,
  tol,
  verbose,
  omega_l,
  log_pi_l,
  omega_bar_l,
  log_omega_bar_l,
  post_params
)
```

## Arguments

- y:

  A numeric matrix. \\T \times d\\ matrix of observations.

- mu_0:

  A numeric vector. Vector of intercept parameters.

- fit_intercept:

  A logical. If `fit_intercept == TRUE`, then an intercept is estimated
  and `mu_0` gets updated.

- refit:

  A logical. If `refit == TRUE`, then the MICH algorithm is initialized
  by the fit provided in `post_params`, otherwise the null model
  \\\boldsymbol{\mu}\_t = \mathbf{0}\\ and \\\bar{\pi}\_{\ell t} = 1/T\\
  is used as the initialization.

- max_iter:

  An integer. Maximum number of iterations. If ELBO does not converge
  before `max_iter` is reached, then `converged == FALSE` in the
  returned fit object.

- tol:

  A scalar. Convergence tolerance for relative increase in ELBO.

- verbose:

  A logical. If `verbose == TRUE`, then the value of the ELBO is printed
  every 5000th iteration.

- omega_l:

  A scalar. Prior precision parameter for mean-scp components of model.

- log_pi_l:

  A numeric matrix. A \\T \times L\\ matrix of prior log change-point
  location probabilities for each of the \\L\\ mean change-points.

- omega_bar_l, log_omega_bar_l:

  Numeric vectors. Vector of posterior precision parameters
  \\\\\bar{\omega}\_t\\\_{t=1}^T\\ and log evaluations such that
  \\V(\mathbf{b}\_\ell\|\tau=t) = \bar{\omega}\_t\mathbf{I}\_d\\.

- post_params:

  A list. A length \\L\\ list of the posterior parameters of each
  mean-scp component. Each element of the list is a list containing a
  \\T\times L\\ matrix of posterior mean parameters and length \\T\\
  vectors of posterior change-point location probabilities and their log
  evaluations.

## Value

A List. Parameters of the variational approximation the MICH posterior
distribution, including:

- `L`: An integer. Number of components included in model.

- `residual`: A numeric matrix. Residual \\\mathbf{r}\_{1:T}\\ after
  subtracting out each \\E\[\boldsymbol{\mu}\_{\ell t}\]\\ from
  \\\mathbf{y}\_{1:T}\\.

- `mu`: A numeric matrix. Posterior estimate of \\\Sigma\_{\ell=1}^L
  E\[\boldsymbol{\mu}\_{\ell,1:T}\|\mathbf{y}\_{1:T}\]\\.

- `mu_0`: A numeric vector. Estimate of the intercept.

- `post_params`: A list. List of posterior parameters for each mean-scp
  component.

- `elbo`: A numeric vector. Value of the ELBO after each iteration.

- `converged`: A logical. Indicates whether relative increase in the
  ELBO is less than `tol`.
