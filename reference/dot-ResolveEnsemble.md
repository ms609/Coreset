# Resolve an expanded ensemble into the winning subset

Shared tail of the two ensemble drivers: solves the specs' distinct
seeds in one batch via the driver's `RunPasses` closure, then returns
the subset maximising \\T_k\\. The returned vector carries the
`strategy_results` (one record per label) and `winning_strategy` (all
tied-best labels) attributes.

## Usage

``` r
.ResolveEnsemble(expanded, labels, RunPasses)
```

## Arguments

- expanded:

  List of `list(label, s1)` specs from
  [`.ExpandAnchors()`](https://ms609.github.io/MaxMin/reference/dot-ExpandAnchors.md).

- labels:

  Character vector of labels (one per spec).

- RunPasses:

  Closure mapping a vector of distinct seeds to a list of
  `list(idx, tK)`, one per seed and in the same order.

## Value

`.ResolveEnsemble()` returns an integer vector of selected indices with
attributes.
