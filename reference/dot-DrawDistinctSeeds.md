# Draw distinct furthest-point seeds from random pivots

Used by
[`.GonzEnsemble()`](https://ms609.github.io/MaxMin/reference/dot-GonzEnsemble.md)
and
[`.GonzEnsembleFromPoints()`](https://ms609.github.io/MaxMin/reference/dot-GonzEnsembleFromPoints.md)
to expand the `"random_furthest"` token (see
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)):
walks distinct random pivots (a partial shuffle of `1:nPts`, so no pivot
is ever tried twice), resolves each pivot's furthest-point seed via
`SeedFunc`, and collects distinct seed indices until `nSeeds` are found
or the draw budget is spent. A `maxDraws` cap bounds the work when the
reachable seed pool is smaller than `nSeeds`. Returns between 1 and
`nSeeds` distinct indices (ascending); set a seed
([`set.seed()`](https://rdrr.io/r/base/Random.html)) for a reproducible
set.

## Usage

``` r
.DrawDistinctSeeds(SeedFunc, nPts, nSeeds, maxDraws = NULL)
```

## Arguments

- SeedFunc:

  Function mapping a pivot index to its furthest-point seed index.

- nPts:

  Integer number of points.

- nSeeds:

  Integer target number of distinct seeds (`>= 1`).

- maxDraws:

  Integer cap on the number of distinct pivots tried. Default
  `max(40 * nSeeds, 100)`.

## Value

`.DrawDistinctSeeds()` returns an integer vector of distinct seed
indices (length in `[1, nSeeds]`).
