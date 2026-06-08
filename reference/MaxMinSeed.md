# Peripheral seed index for Gonzalez farthest-first selection

Returns the index of a single deterministic peripheral seed, the
starting point used by
[`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md)
under the corresponding `seed` strategy. Useful when composing a custom
selection pass.

## Usage

``` r
MaxMinSeed(
  d = NULL,
  points = NULL,
  method = c("peripheral", "centroid", "random_furthest", "diameter", "anti_medoid",
    "medoid", "rowsum", "rownorm")
)
```

## Arguments

- d:

  A `dist` object or square symmetric numeric matrix. Ignored when
  `points` is supplied.

- points:

  Optional `N x dim` numeric coordinate matrix; when supplied the seed
  is computed from coordinates in `O(N)` memory. Required for the
  `"centroid"` anchor, which has no distance-matrix form.

- method:

  Anchor name (see above). Default `"peripheral"`.

## Value

Integer seed index in `[1, N]`.

## Details

Anchors:

- `"centroid"`:

  The point farthest from the coordinate mean (`argmax ||x - x_bar||`),
  an `O(N * dim)` approximate diameter endpoint. Computed from
  coordinates, so it requires `points`; it is unavailable on the
  distance-matrix path, where `"peripheral"` serves the same role.

- `"peripheral"`:

  Two sweeps: the point furthest from point 1, then the point furthest
  from that (a diameter-endpoint approximation), in `O(N)`. The only
  anchor reachable from a distance-column oracle (the function path of
  [`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md)).

- `"random_furthest"`:

  The point furthest from a random pivot, in `O(N)`. The pivot is drawn
  from a fixed internal seed, isolated from the ambient RNG, so the
  index is reproducible across sessions and machines.

- `"diameter"`:

  A row endpoint of the diameter pair (the maximum pairwise distance).
  Degenerate data (`d_max <= 0`) falls back to 1.

- `"anti_medoid"`:

  The point furthest from the 1-median (medoid).

- `"medoid"`:

  The 1-median itself, `which.min(rowSums(d))`.

- `"rowsum"`:

  The point maximising the sum of distances to all others (the
  1-anti-median).

- `"rownorm"`:

  The point maximising `sqrt(sum d^2)`, the L2 counterpart of
  `"rowsum"`.

## See also

[`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md),
which seeds and runs the greedy pass in one call.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
MaxMinSeed(d, method = "diameter")
#> [1] 14
Gonzalez(d, 5L, seed = MaxMinSeed(d, method = "diameter"))
#> [1] 14  4 26  5 28
```
