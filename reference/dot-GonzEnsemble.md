# Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)

Runs Gonzalez from each requested peripheral anchor and returns the
subset maximising \\T_k\\. Internal driver for the ensemble path of
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
(triggered when `strategy` is a character vector of length \> 1 or
`"random_furthest"`). The `"random_furthest"` token draws `nSeeds`
distinct furthest-point seeds via
[`.DrawDistinctSeeds()`](https://ms609.github.io/MaxMin/reference/dot-DrawDistinctSeeds.md).
The returned vector carries `strategy_results` and `winning_strategy`
(character vector of all tied-best strategies, with random starts
labelled `random_furthest1`, `random_furthest2`, ...) attributes.

## Usage

``` r
.GonzEnsemble(d, m, anchors = "peripheral", nSeeds = 3L)
```

## Arguments

- d:

  Square numeric distance matrix (already coerced).

- m:

  Integer subset size (`1 <= m < nrow(d)`).

- anchors:

  Character vector of anchor names.

- nSeeds:

  Integer number of distinct random-furthest seeds to draw when
  `"random_furthest"` is in `anchors`.

## Value

`.GonzEnsemble()` returns an integer vector of selected indices with
attributes.
