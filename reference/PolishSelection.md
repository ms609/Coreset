# Local-search polish for a max-min diversity selection

Given a selection `idx` of size `k >= 2`, search for 1-swaps that
increase \\T_k = \min\_{a \neq b \in S} d(a, b)\\ (or, at fixed `T_k`,
reduce the number of pairs achieving it). The neighbourhood is
critical-edge-anchored: the candidate removal `x` ranges over endpoints
of pairs achieving the current minimum, and the candidate insertion `w`
over the `limit` nearest neighbours of `x` in the full data (excluding
the current selection).

## Usage

``` r
PolishSelection(
  d,
  idx,
  limit = 20L,
  max_passes = 200L,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- d:

  A `dist` object or square symmetric numeric matrix of pairwise
  distances.

- idx:

  Integer vector of selected row/col indices (`k = length(idx)`).

- limit:

  Integer: maximum neighbour rank to scan per critical endpoint. Default
  `20L`. Smaller values are faster but may miss improving swaps.

- max_passes:

  Integer: hard cap on outer iterations (safety guard; theoretical
  termination is guaranteed by the lex-monotone objective). Default
  `200L`.

- progress:

  Logical; show a start/done status line. Default: `TRUE` in interactive
  sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`).

## Value

Integer vector of the same length as `idx`, with attributes `"passes"`
and `"swaps"` (both equal to the number of accepted swaps).

## Details

Acceptance follows Della Croce et al.'s flat-landscape rule: accept a
swap if `new_T > T` *or* `new_T == T` *and* the count of pairs at the
minimum strictly decreases. Both `T_k` and `-n_critical` are monotone
non-decreasing across accepted swaps, so iteration terminates without a
tabu list.

The intended use is as a post-step on any greedy max-min selector — for
example
[`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md).
Calling `PolishSelection()` on a 1-swap local optimum is a no-op.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d   <- dist(pts)
s0  <- Gonzalez(d, 5L, seed = 1L)
s1  <- PolishSelection(d, s0)
TkScore(d, s1) >= TkScore(d, s0)
#> [1] TRUE
```
