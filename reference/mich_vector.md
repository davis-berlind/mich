# Univariate Multiple Independent Change-Point (MICH) Model

Fits the univariate version of the MICH model for mean and variance
change-points. Number of change-points can either be fixed or
`mich_vector()` will search for the number of changes that maximizes the
ELBO when `(L_auto | K_auto | J_auto) == TRUE`.

## Usage

``` r
mich_vector(
  y,
  fit_intercept,
  fit_scale,
  standardize,
  J,
  L,
  K,
  J_auto,
  L_auto,
  K_auto,
  J_max,
  L_max,
  K_max,
  pi_j_weighted,
  pi_l_weighted,
  pi_k_weighted,
  tol,
  verbose,
  max_iter,
  reverse,
  detect,
  merge_level,
  merge_prob,
  restart,
  n_search,
  increment,
  omega_j,
  u_j,
  v_j,
  log_pi_j,
  omega_l,
  log_pi_l,
  u_k,
  v_k,
  log_pi_k
)
```

## Arguments

- y:

  A numeric vector. Length \\T\\ vector of observations.

- fit_intercept:

  A logical. If `fit_intercept == TRUE`, then an intercept is estimated,
  otherwise it is assumed that \\\mu_0 = 0\\.

- fit_scale:

  A logical. If `fit_scale == TRUE`, then the initial precision is
  estimated, otherwise it is assumed that \\\lambda_0 = 1\\.

- standardize:

  A logical. If `standardize == TRUE`, then `y` is centered and rescaled
  before fitting.

- L, K, J:

  Integers. Respective number of mean-scp, var-scp, and meanvar-scp
  components included in the model. If `L_auto == TRUE`,
  `K_auto == TRUE`, or `J_auto == TRUE` then `L`, `K`, and `J` lower
  bound the number of each kind of change-point in the model.

- L_auto, K_auto, J_auto:

  Logicals. If `L_auto == TRUE`, `K_auto == TRUE`, and/or
  `J_auto == TRUE`, then `mich_vector()` returns the \\L\\ between `L`
  and `L_max`, the \\K\\ between `K` and `K_max`, and/or the \\J\\
  between `J` and `J_max` that maximize the ELBO (see Appendix C.4 of
  Berlind, Cappello, Madrid Padilla (2025)).

- L_max, K_max, J_max:

  Integers. If `L_auto == TRUE`, `K_auto == TRUE`, or `J_auto == TRUE`
  then `L_max`, `K_max`, and `J_max` upper bound the number of each kind
  of change-point in the model.

- pi_l_weighted, pi_k_weighted, pi_j_weighted:

  Logicals. If `TRUE`, then the weighted priors specified in Appendix
  C.2 of Berlind, Cappello, Madrid Padilla (2025) are used.

- tol:

  A scalar. Convergence tolerance for relative increase in ELBO.

- verbose:

  A logical. If `verbose == TRUE` and `L_auto == FALSE`,
  `K_auto == FALSE`, and `J_auto == FALSE` then the value of the ELBO is
  printed every 5000th iteration. If `verbose == TRUE` and any of
  `L_auto`, `K_auto`, or `J_auto` are `TRUE`, then the value of the ELBO
  is printed for each combinatin of \\(L,K,J)\\ as `mich_vector()`
  searches for the parameterization that maximized the EBLO.

- max_iter:

  An integer. Maximum number of iterations. If ELBO does not converge
  before `max_iter` is reached, then `converged == FALSE` in the
  returned fit object.

- reverse:

  A logical. If `reverse == TRUE` then MICH is fit to
  \\\mathbf{y}\_{T:1}\\ and the model parameters are reversed in
  post-processing.

- detect:

  A scalar. The detection criteria. The \\i^{\text{th}}\\ component of
  the model detects a change-point only if the posterior credible set
  for that component contains fewer than `detect` indices.

- merge_level:

  A scalar. A value between (0,1) for the significance level to
  construct credible sets at when merging. A model component is only
  considered to be a candidate for merging if its `merge_level`-level
  credible set contains fewer than `detect` indices.

- merge_prob:

  A scalar. A value between (0,1) indicating the merge criterion. If the
  posterior probability that two components identify the same change is
  greater than `merge_level`, then those components are merged (see
  Appendix C.3 of Berlind, Cappello, Madrid Padilla (2025)).

- restart:

  A logical. If `restart == TRUE` and `L_auto`, `K_auto`, or `J_auto`
  are `TRUE`, then after `n_search` increments of `L`, `K`, and/or `J`,
  if the ELBO has not increased, `mich_vector()` will restart by setting
  the `L`, `K`, and `J` components to the null model initialization
  (except for the components with maximum posterior probabilities \>
  0.9) then refit and begin the search again.

- n_search:

  An integer. Grid search parameter. Number of times to increment `L`,
  `K`, and/or `J` before terminating automatic procedure when `L_auto`,
  `K_auto`, or `J_auto` are `TRUE`.

- increment:

  An integer. Number of components to increment `L`, `K` and `J` by when
  `L_auto`, `K_auto`, or `J_auto` are `TRUE`.

- omega_j:

  A scalar. Prior precision parameter for meanvar-scp components of the
  model.

- u_j:

  A scalar. Prior shape parameter for meanvar-scp components of the
  model.

- v_j:

  A scalar. Prior rate parameter for meanvar-scp components of the
  model.

- log_pi_j:

  A numeric matrix. A \\T \times J\\ matrix of prior log change-point
  location probabilities for each of the \\J\\ mean-variance
  change-points.

- omega_l:

  A scalar. Prior precision parameter for mean-scp components of the
  model.

- log_pi_l:

  A numeric matrix. A \\T \times L\\ matrix of prior log change-point
  location probabilities for each of the \\L\\ mean change-points.

- u_k:

  A scalar. Prior shape parameter for var-scp components of the model.

- v_k:

  A scalar. Prior rate parameter for var-scp components of the model.

- log_pi_k:

  A numeric matrix. A \\T \times K\\ matrix of prior log change-point
  location probabilities for each of the \\K\\ variance change-points.

## Value

A list. Parameters of the variational approximation the MICH posterior
distribution, including:

- `y`: A numeric vector. Original data.

- `L`,`K`,`J`: Integers. Number of mean-scp, var-scp, and meanvar-scp
  components included in the model.

- `residual`: A numeric vector. Residual \\\tilde{\mathbf{r}}\_{1:T}\\
  (see (B.4) of Berlind, Cappello, Madrid Padilla (2025)).

- `mu`: A numeric vector. Posterior estimate of \\\Sigma\_{\ell=1}^L
  E\[\mu\_{\ell,1:T}\|\mathbf{y}\_{1:T}\] + \Sigma\_{j=1}^J
  E\[\mu\_{j,1:T}\|\mathbf{y}\_{1:T}\]\\.

- `lambda`: A numeric vector. Posterior estimate of \\\Pi\_{k=1}^K
  E\[\lambda\_{k,1:T}\|\mathbf{y}\_{1:T}\] \times \Pi\_{j=1}^J
  E\[\lambda\_{j,1:T}\|\mathbf{y}\_{1:T}\]\\.

- `delta`: A numeric vector. Posterior estimate of equation (B.4) of
  Berlind, Cappello, Madrid Padilla (2025).

- `mu_0`: A scalar. Estimate of the intercept.

- `lambda_0`: A scalar. Estimate of the initial precision.

- `elbo`: A numeric vector. Value of the ELBO after each iteration.

- `converged`: A logical. Indicates whether relative increase in the
  ELBO is less than `tol`.

- `meanvar_model`: A list. List of meanvar-scp posterior parameters:

  - `pi_bar`: A numeric matrix. A \\T \times J\\ matrix of posterior
    change-point location probabilities.

  - `b_bar`: A numeric matrix. A \\T \times J\\ matrix of posterior mean
    parameters.

  - `omega_bar`: A numeric matrix. A \\T \times J\\ matrix of posterior
    precision parameters.

  - `v_bar`: A numeric matrix. A \\T \times J\\ matrix of posterior rate
    parameters.

  - `u_bar`: A numeric vector. A length \\T\\ vector of posterior shape
    parameters.

  - `mu_lambda_bar`: A numeric matrix. A \\T \times J\\ matrix of scaled
    posterior mean signals
    \\E\[\lambda\_{jt}\mu\_{jt}\|\mathbf{y}\_{1:T}\]\\.

  - `mu2_lambda_bar`: A numeric matrix. A \\T \times J\\ matrix of
    scaled posterior squared mean signals
    \\E\[\lambda\_{jt}\mu^2\_{jt}\|\mathbf{y}\_{1:T}\]\\.

  - `lambda_bar`: A numeric matrix. A \\T \times J\\ matrix of posterior
    precision signals \\E\[\lambda\_{jt}\|\mathbf{y}\_{1:T}\]\\.

- `mean_model`: A list. List of mean-scp posterior parameters:

  - `pi_bar`: A numeric matrix. A \\T \times L\\ matrix of posterior
    change-point location probabilities.

  - `b_bar`: A numeric matrix. A \\T \times L\\ matrix of posterior mean
    parameters.

  - `omega_bar`: A numeric matrix. A \\T \times L\\ matrix of posterior
    precision parameters.

  - `mu_bar`: A numeric matrix. A \\T \times L\\ matrix of posterior
    mean signals \\E\[\mu\_{\ell t}\|\mathbf{y}\_{1:T}\]\\.

  - `mu2_bar`: A numeric matrix. A \\T \times L\\ matrix of posterior
    squared mean signals \\E\[\mu\_{\ell t}^2\|\mathbf{y}\_{1:T}\]\\.

- `var_model`: A list. List of var-scp posterior parameters:

  - `pi_bar`: A numeric matrix. A \\T \times K\\ matrix of posterior
    change-point location probabilities.

  - `v_bar`: A numeric matrix. A \\T \times K\\ matrix of posterior rate
    parameters.

  - `u_bar`: A numeric vector. A length \\T\\ vector of posterior shape
    parameters.

  - `lambda_bar`: A numeric matrix. A \\T \times K\\ matrix of posterior
    precision signals \\E\[\lambda\_{kt}\|\mathbf{y}\_{1:T}\]\\.
