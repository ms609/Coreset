# Resolve an expanded ensemble into the winning subset

Shared tail of the two ensemble drivers: solves each expanded spec via
the driver's `RunGonz` closure (which deduplicates repeated seeds
through its own cache), then returns the subset maximising \\T_k\\. The
returned vector carries the `strategy_results` (one record per label)
and `winning_strategy` (all tied-best labels) attributes.

## Usage

``` r
.ResolveEnsemble(expanded, labels, RunGonz)
```

## Arguments

- expanded:

  List of `list(label, s1)` specs from
  [`.ExpandAnchors()`](https://ms609.github.io/MaxMin/reference/dot-ExpandAnchors.md).

- labels:

  Character vector of labels (one per spec).

- RunGonz:

  Closure mapping a seed `s1` to `list(idx, tK)`.

## Value

Integer vector of selected indices with attributes.
