# Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)

Runs Gonzalez from each requested deterministic peripheral anchor and
returns the subset maximising \\T_k\\. Internal driver for the ensemble
path of
[`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md)
(triggered when `seed` is a character vector of length \> 1). The
returned vector carries `strategy_results` and `winning_strategy`
(character vector of all tied-best strategies) attributes.

## Usage

``` r
.GonzEnsemble(
  d,
  n,
  anchors = c("diameter", "anti_medoid", "rowsum", "rownorm")
)
```

## Arguments

- d:

  Square numeric distance matrix (already coerced).

- n:

  Integer subset size (`1 <= n < nrow(d)`).

- anchors:

  Character vector of anchor names.

## Value

Integer vector of selected indices with attributes.
