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
  on the selected sub-coordinates only (`k x k`).

## Value

`MinDist()` returns a numeric specifying the minimum distance between
two selected points (or `NA_real_`, if `length(idx) < 2`).

## Details

The solvers in this package
([`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)) already
attach the achieved \\T_k\\ as a `score` attribute. `MinDist()` allows
arbitrary selections to be scored, or an existing selection to be scored
against a different distance matrix.

## See also

[`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md) and
[`ExactMaxMin()`](https://ms609.github.io/Coreset/reference/ExactMaxMin.md).

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
MinDist(d, FarFirst(5L, d))
#> [1] 1.765223
```
