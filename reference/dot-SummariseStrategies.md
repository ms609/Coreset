# Print the per-strategy \\T_k\\ table of a [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) ensemble

One row per strategy tried, ordered best (largest \\T_k\\) first, with
each tied-best strategy marked `*`. A bare single pass (no
`strategy_results`) produces nothing.

## Usage

``` r
.SummariseStrategies(object)
```

## Arguments

- object:

  A `MaxMinSelection` from an ensemble
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
  call.

## Value

Invisibly `NULL`; called for the side effect.
