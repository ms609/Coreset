# Peripheral seed index for Gonzalez selection (coordinates)

Coordinate counterpart of
[`.PickPoint()`](https://ms609.github.io/MaxMin/reference/dot-PickPoint.md);
each anchor is computed from the `…FromPoints_cpp` primitives,
bit-identical to the matrix path on Euclidean data.

## Usage

``` r
.PickPoints(points, strategy)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- strategy:

  Anchor name; see
  [`PickPoint()`](https://ms609.github.io/MaxMin/reference/PickPoint.md).
  Also accepts `"first"` (1).

## Value

`.PickPoints()` returns an integer seed index.
