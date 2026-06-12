# Minimum pairwise distance within a selection, from coordinates

Coordinate counterpart of `.SubsetScore(d, idx, "min_pairwise")`.
Computes [`stats::dist()`](https://rdrr.io/r/stats/dist.html) on the
selected sub-coordinates only (`k x k`, never the full matrix); the
per-pair bits are identical to the corresponding entries of the full
distance matrix, so the returned scalar matches the matrix path.

## Usage

``` r
.MinPairwiseFromPoints(points, idx)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- idx:

  Integer indices of the selection.

## Value

`.MinPairwiseFromPoints()` returns a numeric scalar; `NA_real_` if
`length(idx) < 2`.
