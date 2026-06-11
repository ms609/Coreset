# Minimum pairwise distance within a selection

Returns the minimum pairwise distance among selected points \\T_k\\. A
set of points that is more dispersed will exhibit a higher value.

## Usage

``` r
MinDist(d = NULL, idx, points = NULL)
```

## Arguments

- d:

  Pairwise distance matrix or `dist` object. Ignored when `points` is
  supplied.

- idx:

  Integer vector of selected row/col indices.

- points:

  Optional `N x dim` numeric coordinate matrix. When supplied, the score
  is computed from [`stats::dist()`](https://rdrr.io/r/stats/dist.html)
  on the selected sub-coordinates only (`k x k`), never the full `N x N`
  matrix (`d` is then unused). For Euclidean data the result is
  identical to the matrix path.

## Value

Numeric scalar; `NA_real_` if `length(idx) < 2`.

## Details

The solvers in this package
([`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md)) already
attach the achieved \\T_k\\ as a `score` attribute, so `MinDist()` is
mainly for scoring a selection produced elsewhere – a matrix-free or
externally generated index set – or for re-scoring an existing selection
against a different distance matrix.

## See also

[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) and
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md),
whose results already carry the objective.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
MinDist(d, FarFirst(d, 5L))
#> [1] 1.765223
```
