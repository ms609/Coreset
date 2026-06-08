# Minimum pairwise distance within a selection (T_k = k-centre objective)

Returns the minimum pairwise distance among selected points, the
canonical k-centre objective \\T_k\\. Higher values indicate a more
spread-out (better) selection.

## Usage

``` r
TkScore(d = NULL, idx, points = NULL)
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

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
TkScore(d, Gonzalez(d, 5L))
#> [1] 1.765223
```
