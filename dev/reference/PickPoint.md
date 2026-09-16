# Seed to initialize farthest-first selection

`PickPoint()` implements a range of strategies to select a seed for
greedy farthest-first selection. Propitious seeds yield better
solutions.

## Usage

``` r
PickPoint(
  d = NULL,
  points = NULL,
  strategy = c("peripheral", "anti_centroid", "random_furthest", "diameter",
    "anti_medoid", "medoid", "rowsum", "rownorm"),
  N = NULL,
  nSeeds = 1L
)
```

## Arguments

- d:

  A `dist` object, a square symmetric numeric matrix, or a distance
  function as accepted by
  [`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md).
  Ignored when `points` is supplied.

- points:

  Optional `N x dim` numeric coordinate matrix; when supplied the seed
  is computed from coordinates in `O(N)` memory. Required for the
  `"anti_centroid"` anchor, which has no distance-matrix form.

- strategy:

  Character specifying method to employ:

  `"peripheral"`(default)

  :   Two sweeps: the point furthest from point 1, then the point
      furthest from that (a diameter-endpoint approximation). \\O(N)\\.

  `"anti_centroid"`

  :   The point farthest from the coordinate mean (\\\arg\max \\x -
      \bar{x}\\\\). \\O(N \* dim)\\. Requires `points`.

  `"random_furthest"`

  :   The point furthest from a random pivot. \\O(N)\\ per pivot.

  `"diameter"`

  :   A row endpoint of the diameter pair (the maximum pairwise
      distance).

  `"medoid"`

  :   The 1-median (medoid): the point minimising the sum of distances
      to all others.

  `"anti_medoid"`

  :   The point furthest from the 1-median (medoid).

  `"rowsum"`

  :   The point maximising the sum of distances to all others (the
      1-anti-median).

  `"rownorm"`

  :   The point maximising \\\sqrt(\sum{d^2})\\, the L2 counterpart of
      `"rowsum"`.

  Only `"peripheral"` and `"random_furthest"` are supported when `d` is
  a function.

- N:

  Integer: the total number of elements. Required (and used) only on the
  distance-column oracle path, where it cannot be inferred from the
  closure; ignored for the matrix and coordinate paths.

- nSeeds:

  Integer specifying how many distinct seeds to obtain under
  `"random_furthest"`; more pivots are evaluated until `nSeeds` distinct
  seeds are found.

## Value

`PickPoint()` returns an integer that identifies the index of a proposed
seed in `d` or `points`; under `"random_furthest"` with `nSeeds > 1`, up
to `nSeeds` distinct indices, in ascending order.

## See also

[`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md),
which seeds and runs the greedy pass in one call.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
PickPoint(d, strategy = "diameter")
#> [1] 14
FarFirst(5L, d, strategy = PickPoint(d, strategy = "diameter"))
#> 5 elements (14 4 26 5 28) selected by farthest-first, each at distance >= 1.765

# Seeds for three starts, from distances computed one column at a time.
# The seeds drawn are identical to those selected by FarFirst().
Column <- function(i) sqrt(colSums((t(pts) - pts[i, ]) ^ 2))
set.seed(2)
PickPoint(Column, strategy = "random_furthest", N = nrow(pts), nSeeds = 3)
#> [1] 14 24 26
set.seed(2)
FarFirst(5L, Column, N = nrow(pts), strategy = "random_furthest", nSeeds = 3)
#> 5 elements (14 4 26 5 28) selected by farthest-first (best of 3 strategies, winner random_furthest1), each at distance >= 1.765
```
