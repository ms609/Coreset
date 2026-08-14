# Mean dispersion of a selection

`MeanDist()` reports the sum of pairwise distances divided by the number
of selected elements, \$\$f(S) = \frac{\displaystyle\sum\_{i \< j,\\ i,j
\in S} d\_{ij}}{\|S\|}\$\$, the objective maximised by
[`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md).

## Usage

``` r
MeanDist(d, idx)
```

## Arguments

- d:

  Pairwise distance matrix or `dist` object.

- idx:

  Integer vector of selected row/col indices.

## Value

`MeanDist()` returns a numeric scalar, or `NA_real_` if
`length(idx) < 2`.

## See also

[`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md) which
maximises this objective;
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) for
the max-min (MMDP) analogue.

## Examples

``` r
# The max-mean problem is defined for signed dissimilarities; with these the
# optimal subset has an interior size, chosen to maximise mean dispersion.
set.seed(1)
x <- matrix(runif(100, -5, 5), 10)
d <- (x + t(x)) / 2          # symmetric, signed
selection <- MaxMean(d)
selection                    # 5 of the 10 elements: {2, 5, 6, 8, 10}
#> 5 elements (2 5 6 8 10) selected by MaxMean RLTS, f = 2.714
MeanDist(d, selection)       # 2.714 — equals attr(selection, "score")
#> [1] 2.714473
```
