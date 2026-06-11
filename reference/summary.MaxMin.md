# Multi-line summaries of MaxMin solver results

A fuller counterpart to
[print.MaxMinSelection](https://ms609.github.io/MaxMin/reference/print.MaxMin.md):
the one-line headline, followed by the achieved objective(s), search
effort, and – for a
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
ensemble – the per-strategy \\T_k\\ table.

## Usage

``` r
# S3 method for class 'MaxMinSelection'
summary(object, ...)

# S3 method for class 'MaxMinExact'
summary(object, ...)
```

## Arguments

- object:

  A `MaxMinSelection` or `MaxMinExact` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

`object`, invisibly.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
summary(FarFirst(dist(pts), 5L))
#> 5 elements (14 4 26 5 28) selected by Gonzalez farthest-first (best of 3 strategies, winner random_furthest2), each at distance >= 1.765
#>   strategies tried (3), best marked *:
#>     strategy           seed  T_k
#>   * random_furthest2     14  1.765
#>     random_furthest1     24  1.701
#>     random_furthest3     26  1.468
```
