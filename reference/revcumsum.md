# Reverse cumulative sum

`revcumsum()` takes a length \\n\\ vector \\\\x_i\\\_{i=1}^n\\ and
returns \\\\\Sigma\_{j=i}^n x_j\\\_{i=1}^n\\.

## Usage

``` r
revcumsum(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A numeric vector. Reverse cumulative sum of the elements of `x`.

## Examples

``` r
revcumsum(1:10)
#>  [1] 55 54 52 49 45 40 34 27 19 10
```
