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
```

## Arguments

- object:

  A `MaxMinSelection` object (from any solver).

- ...:

  Ignored; present for S3 compatibility.

## Value

`summary.MaxMin()` returns `object`, invisibly.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`print.MaxEntropy`](https://ms609.github.io/MaxMin/reference/print.MaxEntropy.md),
[`print.MaxMeanSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMeanSelection.md),
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md),
[`print.MaxSum`](https://ms609.github.io/MaxMin/reference/print.MaxSum.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
summary(FarFirst(5L, dist(pts)))
#> 5 elements (14 4 26 5 28) selected by farthest-first (best of 3 strategies, winner random_furthest1), each at distance >= 1.765
#>   strategies tried (3), best marked *:
#>     strategy           seed  T_k
#>   * random_furthest1     14  1.765
#>     random_furthest2     24  1.701
#>     random_furthest3     26  1.468
```
