# Multivariate Multiple Independent Change-Point (MICH) Model

Fits the multivariate version of the MICH model for mean change-points.
Number of change-points can either be fixed or `mich_matrix()` will
search for the number of changes that maximizes the ELBO when
`L_auto == TRUE`.

## Usage

``` r
mich_matrix(
  y,
  fit_intercept,
  fit_scale,
  standardize,
  L,
  L_auto,
  L_max,
  pi_l_weighted,
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
  omega_l,
  log_pi_l
)
```

## Arguments

- y:

  A numeric matrix. T x d matrix of observations.

- fit_intercept:

  A logical. If `fit_intercept == TRUE`, then an intercept is estimated,
  otherwise `mu_0 = rep(0,d)`.

- fit_scale:

  A logical. If `fit_scale == TRUE`, then the precision matrix is
  estimated using the inverse of `var(diff(y))`, otherwise it is assumed
  that the precision matrix is equal to `diag(d)`.

- standardize:

  A logical. If `standardize == TRUE`, then `y` is centered and rescaled
  before fitting.

- L:

  An integer. Number of mean-scp components included in model. If
  `L_auto == TRUE` then `L` lower bounds the number of change-points in
  the model.

- L_auto:

  A logical. If `L_auto == TRUE`, then `mich_matrix()` returns the
  number of changes between `L` and `L_max` that maximizes the ELBO (see
  Appendix C.4 of Berlind, Cappello, Madrid Padilla (2025)).

- L_max:

  L An integer. If `L_auto == TRUE` then `L_max` upper bounds the number
  of change-points included in the model.

- pi_l_weighted:

  A logical. If `pi_l_weighted == TRUE`, then the weighted priors
  specified in Appendix C.2 of Berlind, Cappello, Madrid Padilla (2025)
  are used.

- tol:

  A scalar. Convergence tolerance for relative increase in ELBO.

- verbose:

  A logical. If `verbose == TRUE` and `L_auto == FALSE`, then the value
  of the ELBO is printed every 5000th iteration. If `verbose == TRUE`
  and `L_auto == TRUE`, the the value of the ELBO is printed for each L
  as `mich_matrix()` searches over \[`L`, `L_max`\].

- max_iter:

  An integer. Maximum number of iterations. If ELBO does not converge
  before `max_iter` is reached, then `converged == FALSE` in the
  returned fit object.

- reverse:

  A logical. If `reverse == TRUE` then MICH is fit to `y[T:1,]` and the
  model parameters are reversed in post-processing.

- detect:

  A scalar. The detection criteria. Each component of the model detects
  a change-point only if the posterior credible set for that component
  contains fewer than `detect` indices.

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

  A logical. If `restart == TRUE` and `L_auto == TRUE` then after
  `n_search` increments of `L`, if the ELBO has not increased,
  `mich_matrix()` will restart by setting the `L` components to the null
  model initialization (except for the components with maximum posterior
  probabilities \> 0.9) then refit and begin the search again.

- n_search:

  An integer. Grid search parameter. Number of times to increment `L`
  before terminating automatic procedure when `L_auto == TRUE`.

- increment:

  An integer. Number of components to increment `L` by when
  `L_auto == TRUE`.

- omega_l:

  A scalar. Prior precision parameter for mean-scp components of model.

- log_pi_l:

  A numeric matrix. T x L matrix of prior log change-point location
  probabilities for each of the L mean change-points.

## Value

A list. Parameters of the variational approximation the MICH posterior
distribution, including:

- `y`: A numeric matrix. Original data.

- `Sigma`: A numeric matrix. Estimate of the precision if
  `fit_scale == TRUE`.

- `L`: An integer. Number of components included in model.

- `pi_bar`: A numeric matrix. A T x L matrix of posterior change-point
  location probabilites.

- `residual`: A numeric matrix. Residual `y - mu`.

- `mu`: A numeric matrix. Posterior estimate of mean signal.

- `mu_0`: A numeric vector. Estimate of the intercept.

- `post_params`: A list. List of posterior parameters for each mean-scp
  component.

- `elbo`: A numeric vector. Value of the ELBO after each iteration.

- `converged`: A logical. Indicates whether relative increase in the
  ELBO is less than `tol`.
