# Coordinate (matrix-free) multi-anchor Gonzalez ensemble

Coordinate counterpart of
[`.GonzEnsemble()`](https://ms609.github.io/MaxMin/reference/dot-GonzEnsemble.md);
each anchor seed and the greedy expansion are computed from `points` via
the coordinate primitives, so the returned indices and attributes match
the matrix path on Euclidean data. `pivots` indexes points directly, so
the `"random_furthest"` starts also match the matrix path.

## Usage

``` r
.GonzEnsembleFromPoints(
  points,
  m,
  anchors = .kDefaultEnsemble,
  pivots = integer(0)
)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- m:

  Integer subset size.

- anchors:

  Character vector of anchor names.

- pivots:

  Integer vector of pivot indices the `"random_furthest"` token expands
  over (empty contributes none).

## Value

Integer vector of selected indices with attributes.
