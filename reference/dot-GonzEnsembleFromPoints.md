# Coordinate (matrix-free) multi-anchor Gonzalez ensemble

Coordinate counterpart of
[`.GonzEnsemble()`](https://ms609.github.io/Coreset/reference/dot-GonzEnsemble.md);
each anchor seed and the greedy expansion are computed from `points` via
the coordinate primitives, so the returned indices and attributes match
the matrix path on Euclidean data.

## Usage

``` r
.GonzEnsembleFromPoints(points, m, anchors = "random_furthest", nSeeds = 3L)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- m:

  Integer subset size.

- anchors:

  Character vector of anchor names.

- nSeeds:

  Integer number of distinct random-furthest seeds to draw when
  `"random_furthest"` is in `anchors`.

## Value

`.GonzEnsembleFromPoints()` returns an integer vector of selected
indices with attributes.
