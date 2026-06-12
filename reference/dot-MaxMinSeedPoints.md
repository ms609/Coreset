# Peripheral seed index for Gonzalez selection (coordinates)

Coordinate counterpart of
[`.MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/dot-MaxMinSeed.md);
each anchor is computed from the `…FromPoints_cpp` primitives,
bit-identical to the matrix path on Euclidean data.

## Usage

``` r
.MaxMinSeedPoints(points, strategy)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

- strategy:

  Anchor name; see
  [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md).
  Also accepts `"first"` (1).

## Value

Integer seed index.
