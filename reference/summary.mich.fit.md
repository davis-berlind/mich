# Summary method for mich.fit objects

Prints a summary of the resulting fit from calling
[`mich()`](https://davis-berlind.github.io/mich/reference/mich.md),
including the ELBO for the fitted model, estimated change-points with
`level`-level credible sets.

## Usage

``` r
# S3 method for class 'mich.fit'
summary(object, level = 0.95, max_length = NULL, ...)
```

## Arguments

- object:

  A mich.fit object. Output of running
  [`mich()`](https://davis-berlind.github.io/mich/reference/mich.md) on
  a numeric vector or matrix.

- level:

  A scalar. A number between (0,1) indicating the significance level to
  construct credible sets at.

- max_length:

  An integer. Detection threshold, see
  [`mich_sets()`](https://davis-berlind.github.io/mich/reference/mich_sets.md).
  Equal to `log(T)^2` by default.

- ...:

  Additional arguments to be passed to methods.

## Value

A list. A list of summary quantities including:

- `elbo`: The value of the ELBO for the model.

- `converged`: Indicator for whether the model has converged.

- `level`: The significance level used to construct credible sets.

- `L`,`K`,`J`: The number of mean-scp, var-scp, and meanvar-scp
  components included in the model.

- `mean_cp`,`var_cp`,`meanvar_cp`: Lists with `cp`, the estimated
  change-points, and `sets`, their corresponding `level`-level credible
  sets.
