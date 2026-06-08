# Coerce coordinate input for the on-the-fly (matrix-free) samplers

The coordinate paths require a complete numeric `N x dim` matrix with
`double` storage; the C++ kernels reproduce
[`stats::dist()`](https://rdrr.io/r/stats/dist.html)'s exact Euclidean
bits, which is only defined for complete data.

## Usage

``` r
.AsPointsMatrix(points)
```

## Arguments

- points:

  A numeric matrix (or coercible) of point coordinates.

## Value

A `double` numeric matrix.
