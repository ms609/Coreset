# Gonzalez maximin from coordinates (matrix-free)

Coordinate counterpart of
[`.MaximinFrom()`](https://ms609.github.io/MaxMin/reference/dot-MaximinFrom.md):
greedy furthest-point selection that recomputes each needed distance
column from `points` on the fly, never materialising the `N x N` matrix.
Bit-identical selection to the matrix path on Euclidean data.

## Usage

``` r
.MaximinFromPoints(points, k, first, mask = 0L)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- k:

  Integer: target subsample size (`>= 1`).

- first:

  Integer: index of the first selected point.

- mask:

  Integer 1-based index of a point to forbid from selection (`0L` =
  none); used by the anti-medoid path to exclude the medoid.

## Value

`.MaximinFromPoints()` returns an integer vector of selected indices.
