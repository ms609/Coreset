# Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)

Runs Gonzalez from each requested peripheral anchor and returns the
subset maximising \\T_k\\. Internal driver for the ensemble path of
[`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md)
(triggered when `seed` is a character vector of length \> 1). The
`"random_furthest"` token expands to `n_random` fixed-seed
random-furthest starts. The returned vector carries `strategy_results`
and `winning_strategy` (character vector of all tied-best strategies,
with random starts labelled `random_furthest1`, `random_furthest2`, ...)
attributes.

## Usage

``` r
.GonzEnsemble(d, n, anchors = "peripheral", n_random = 0L)
```

## Arguments

- d:

  Square numeric distance matrix (already coerced).

- n:

  Integer subset size (`1 <= n < nrow(d)`).

- anchors:

  Character vector of anchor names.

- n_random:

  Integer; number of starts the `"random_furthest"` token expands to
  (`0` contributes none).

## Value

Integer vector of selected indices with attributes.
