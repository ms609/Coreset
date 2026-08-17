# Format and print Coreset solver results

Terse summaries of the objects returned by the Coreset solvers.

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
  [`ExactMaxMin()`](https://ms609.github.io/Coreset/reference/ExactMaxMin.md)).

- ...:

  Ignored; present for S3 compatibility.

## Value

`print.Coreset()` returns `x`, invisibly. It is called for its
side-effect of printing `format(x)` to the console. `format.Coreset()`
returns a character string reporting the selection size, the selected
indices, the algorithm (and if applicable strategy or proof status), and
the achieved \\T_k\\.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/Coreset/reference/print.KCentre.md),
[`print.MaxEntropy`](https://ms609.github.io/Coreset/reference/print.MaxEntropy.md),
[`print.MaxMeanSelection()`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md),
[`print.MaxSum`](https://ms609.github.io/Coreset/reference/print.MaxSum.md),
[`summary.Coreset`](https://ms609.github.io/Coreset/reference/summary.Coreset.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
print(FarFirst(5L, dist(pts)))
#> 5 elements (14 4 26 5 28) selected by farthest-first (best of 3 strategies, winner random_furthest1), each at distance >= 1.765
```
