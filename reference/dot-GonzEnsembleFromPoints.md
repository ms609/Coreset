# Coordinate (matrix-free) multi-anchor Gonzalez ensemble

Coordinate counterpart of
[`.GonzEnsemble()`](https://ms609.github.io/MaxMin/reference/dot-GonzEnsemble.md);
each anchor seed and the greedy expansion are computed from `points` via
the coordinate primitives, so the returned indices and attributes match
the matrix path on Euclidean data. The random pivots depend only on `N`
and the fixed seed, so the `"random_furthest"` starts also match the
matrix path.

## Usage

``` r
.GonzEnsembleFromPoints(points, n, anchors = .kDefaultEnsemble, n_random = 0L)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- n:

  Integer subset size.

- anchors:

  Character vector of anchor names.

- n_random:

  Integer; number of starts the `"random_furthest"` token expands to
  (`0` contributes none).

## Value

Integer vector of selected indices with attributes.
