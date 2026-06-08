# Coordinate (matrix-free) multi-anchor Gonzalez ensemble

Coordinate counterpart of
[`.GonzEnsemble()`](https://ms609.github.io/MaxMin/reference/dot-GonzEnsemble.md);
each anchor seed and the greedy expansion are computed from `points` via
the coordinate primitives, so the returned indices and attributes match
the matrix path on Euclidean data.

## Usage

``` r
.GonzEnsembleFromPoints(
  points,
  n,
  anchors = c("diameter", "anti_medoid", "rowsum", "rownorm")
)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- n:

  Integer subset size.

- anchors:

  Character vector of anchor names.

## Value

Integer vector of selected indices with attributes.
