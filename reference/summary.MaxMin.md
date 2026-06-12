# Detailed summaries of MaxMin solver results

A fuller counterpart to
[`print.MaxMinSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md):
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

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
summary(FarFirst(dist(pts), 5L))
#> 5 elements (14 4 26 5 28) selected by Gonzalez farthest-first (best of 5 strategies, 3 tied: random_furthest2, random_furthest4, random_furthest5), each at distance >= 1.765
#>   strategies tried (5), best marked *:
#>     strategy           seed  T_k
#>   * random_furthest2     14  1.765
#>   * random_furthest4      5  1.765
#>   * random_furthest5      4  1.765
#>     random_furthest1     24  1.701
#>     random_furthest3     26  1.468
```
