# Draw distinct furthest-point seeds from random pivots

The `nseeds` counterpart to an explicit `pivots` list (see
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)):
repeatedly draws a random pivot with the session RNG, resolves its
furthest-point seed via `seedFn`, and collects distinct seed indices
until `nseeds` are found. Two bounds stop the loop when the reachable
seed pool is smaller than `nseeds`: a consecutive-miss limit (the pool
is likely exhausted once many draws in a row yield only already-seen
seeds) and an absolute draw budget. Returns between 1 and `nseeds`
distinct indices; set a seed
([`set.seed()`](https://rdrr.io/r/base/Random.html)) for a reproducible
set.

## Usage

``` r
.DrawDistinctSeeds(seedFn, nPts, nseeds, maxDraws = NULL, missLimit = NULL)
```

## Arguments

- seedFn:

  Function mapping a pivot index to its furthest-point seed index.

- nPts:

  Integer number of points.

- nseeds:

  Integer target number of distinct seeds (`>= 1`).

- maxDraws:

  Integer absolute draw budget. Default `max(40 * nseeds, 100)`.

- missLimit:

  Integer consecutive-miss limit. Default `max(8 * nseeds, 30)`.

## Value

Integer vector of distinct seed indices (length in `[1, nseeds]`).
