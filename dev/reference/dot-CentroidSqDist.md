# Squared distance of every point to the coordinate anti_centroid

The `O(N * dim)` basis of the `"anti_centroid"` seed: its argmax is the
point farthest from the coordinate mean, an approximate diameter
endpoint.

## Usage

``` r
.CentroidSqDist(points)
```

## Arguments

- points:

  A `double` `N x dim` coordinate matrix.

## Value

`.CentroidSqDist()` returns a numeric vector of length `N` of squared
distances to the mean.
