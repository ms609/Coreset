# Format and print MaxMin solver results

Terse summaries of the objects returned by the MaxMin solvers.

## Usage

``` r
# S3 method for class 'MaxMinSelection'
format(x, ...)

# S3 method for class 'MaxMinSelection'
print(x, ...)
```

## Arguments

- x:

  A `MaxMinSelection` object (from any solver, including
  [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)).

- ...:

  Ignored; present for S3 compatibility.

## Value

`print.MaxMin()` returns `x`, invisibly. It is called for its
side-effect of printing `format(x)` to the console. `format.MaxMin()`
returns a character string reporting the selection size, the selected
indices, the algorithm (and if applicable strategy or proof status), and
the achieved \\T_k\\.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
print(FarFirst(5L, dist(pts)))
#> 5 elements (4 14 26 5 28) selected by farthest-first (best of 5 strategies, 3 tied: random_furthest1, random_furthest2, random_furthest3), each at distance >= 1.765
```
