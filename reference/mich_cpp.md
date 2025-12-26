# MICH Algorithm

Implementation of Algorithms 1 & 2 in Berlind, Cappello, and Madrid
Padilla (2025). This algorithm takes a sequence of \\T\\ observations
\\\mathbf{y}\_{1:T}\\, and iteratively fits \\L\\ mean-scp models, \\K\\
var-scp models, and \\J\\ meanvar-scp models resulting in a variational
approximation to the posterior distribution of the \\L\\ mean, \\K\\
variance, and \\J\\ mean and variance change-points. The algorithm
terminates once the percentage increase in the ELBO falls below `tol`.

## Usage

``` r
mich_cpp(
  y,
  J,
  L,
  K,
  mu_0,
  lambda_0,
  fit_intercept,
  fit_scale,
  refit,
  max_iter,
  verbose,
  tol,
  omega_j,
  u_j,
  v_j,
  log_pi_j,
  pi_bar_j,
  log_pi_bar_j,
  b_bar_j,
  omega_bar_j,
  u_bar_j,
  v_bar_j,
  lgamma_u_bar_j,
  digamma_u_bar_j,
  omega_l,
  log_pi_l,
  pi_bar_l,
  log_pi_bar_l,
  b_bar_l,
  omega_bar_l,
  u_k,
  v_k,
  log_pi_k,
  pi_bar_k,
  log_pi_bar_k,
  u_bar_k,
  v_bar_k,
  lgamma_u_bar_k,
  digamma_u_bar_k
)
```

## Arguments

- y:

  A numeric vector. Length \\T\\ vector of observations.

- L, K, J:

  Integers. The number of mean, variance, and mean-variance
  change-points to include in the model.

- mu_0:

  A scalar. Intercept parameter initialization.

- lambda_0:

  A scalar. Baseline scale parameter initialization.

- fit_intercept:

  A logical. If `fit_intercept == TRUE`, then an intercept is estimated
  and `mu_0` gets updated.

- fit_scale:

  A logical. If `fit_scale == TRUE`, then an initial scale is estimated
  and `lambda_0` gets updated.

- refit:

  A logical. If `refit == TRUE`, then the MICH algorithm is initialized
  using the provided posterior parameters, otherwise the null model
  \\\mu_t = 0\\, \\\lambda_t = 1\\ and \\\bar{\pi}\_{\ell t} =
  \bar{\pi}\_{k t} = \bar{\pi}\_{j t} = 1/T\\ is used as the
  initialization.

- max_iter:

  An integer. Maximum number of iterations. If ELBO does not converge
  before `max_iter` is reached, then `converged == FALSE` in the
  returned fit object.

- verbose:

  A logical. If `verbose == TRUE`, then value of the ELBO is printed
  every 5000th iteration.

- tol:

  A scalar. Convergence tolerance for relative increase in ELBO.

- omega_j, u_j, v_j:

  Scalars. Prior precision, shape, and rate parameters for meanvar-scp
  components of model.

- log_pi_j:

  A numeric matrix. A \\T \times J\\ matrix of prior log change-point
  location probabilities for each of the \\J\\ mean-variance
  change-points.

- pi_bar_j, log_pi_bar_j:

  Numeric matrices. \\T \times J\\ matrices of initialized posterior
  change-point location probabilities and their log evaluations for the
  \\J\\ mean-variance change-points.

- b_bar_j:

  A numeric matrix. A \\T \times J\\ matrix of initialized posterior
  mean parameters of the \\J\\ mean-variance change-points.

- omega_bar_j:

  A numeric matrix. A \\T\times J\\ matrix of initialized posterior
  precision parameters of the \\J\\ mean-variance change-points.

- u_bar_j, lgamma_u_bar_j, digamma_u_bar_j:

  Numeric vectors. Length \\T\\ vectors of posterior shape parameters
  for the meanvar-scp model components and their log-gamma and digamma
  evaluations.

- v_bar_j:

  A numeric matrix. A \\T \times J\\ matrix of initialized posterior
  rate parameters for the \\J\\ mean-variance change-points.

- omega_l:

  A scalar. Prior precision parameter for mean-scp components of model.

- log_pi_l:

  A numeric matrix. A \\T \times L\\ matrix of prior log change-point
  location probabilities for each of the \\L\\ mean change-points.

- pi_bar_l, log_pi_bar_l:

  Numeric matrices. \\T \times L\\ matrices of initialized posterior
  change-point location probabilities and their log evaluations for the
  \\L\\ mean change-points.

- b_bar_l:

  A numeric matrix. A \\T \times L\\ matrix of initialized posterior
  mean parameters of the \\L\\ mean change-points.

- omega_bar_l:

  A numeric matrix. A \\T\times L\\ matrix of initialized posterior
  precision parameters of the \\L\\ mean change-points.

- u_k, v_k:

  Scalars. Prior shape and rate parameters for var-scp components of
  model.

- log_pi_k:

  A numeric matrix. A \\T \times K\\ matrix of prior log change-point
  location probabilities for each of the \\K\\ variance change-points.

- pi_bar_k, log_pi_bar_k:

  Numeric matrices. \\T \times K\\ matrices of initialized posterior
  change-point location probabilities and their log evaluations for the
  \\K\\ variance change-points.

- u_bar_k, lgamma_u_bar_k, digamma_u_bar_k:

  Numeric vectors. Length \\T\\ vectors of posterior shape parameters
  for the var-scp model components and their log-gamma and digamma
  evaluations.

- v_bar_k:

  A numeric matrix. A \\T \times K\\ matrix of initialized posterior
  rate parameters for the \\K\\ variance change-points.

## Value

A List. Parameters of the variational approximation the MICH posterior
distribution, including:

- `y`: A numeric vector. Original data set.

- `residual`: A numeric vector. Residual \\\mathbf{r}\_{1:T}\\ after
  subtracting out each \\E\[\mu\_{\ell t}\]\\ and \\E\[\lambda\_{j t}
  \mu\_{j t}\]/E\[\lambda\_{j t}\]\\ from \\\mathbf{y}\_{1:T}\\.

- `mu`: A numeric vector. Posterior estimate of \\\Sigma\_{\ell=1}^L
  E\[\mu\_{\ell,1:T}\|\mathbf{y}\_{1:T}\] + \Sigma\_{j=1}^J
  E\[\mu\_{j,1:T}\|\mathbf{y}\_{1:T}\]\\.

- `lambda`: A numeric vector. Posterior estimate of \\\Pi\_{\\k=1}^K
  E\[\lambda\_{k,1:T}\|\mathbf{y}\_{1:T}\] \times \Pi\_{j=1}^J
  E\[\lambda\_{j,1:T}\|\mathbf{y}\_{1:T}\]\\.

- `delta`: A numeric vector. Posterior variance correction term (see
  (B.6) or Berlind, Cappello, and Madrid Padilla (2025)).

- `converged`: A logical. Indicates whether relative increase in the
  ELBO is less than `tol`.

- `elbo`: A numeric vector. Value of the ELBO after each iteration.

- `mu_0`: A scalar. Estimate of the intercept.

- `lambda_0`: A scalar. Estimate of the initial precision.

- `J`, `L`, `K`: Integers. number of mean, variance, and mean-variance
  components.

- `J_model`: A list. A list of the posterior parameters for each of the
  \\J\\ meanvar-scp components (only included in `J > 0`).

- `L_model`: A list. A list of the posterior parameters for each of the
  \\L\\ mean-scp components (only included in `L > 0`).

- `K_model`: A list. A list of the posterior parameters for each of the
  \\K\\ var-scp components (only included in `K > 0`).
