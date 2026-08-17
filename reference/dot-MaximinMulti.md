# Gonzalez maximin from several starting indices at once

Ensemble counterpart of
[`.MaximinFrom()`](https://ms609.github.io/Coreset/reference/dot-MaximinFrom.md)
/
[`.MaximinFromPoints()`](https://ms609.github.io/Coreset/reference/dot-MaximinFromPoints.md):
solves one greedy pass per seed. Each pass is an independent function of
its seed, so under `mc.cores > 1` the passes run concurrently, one per
thread; a lone seed, or a problem large enough for a single pass to
occupy every thread itself, falls back to the per-pass parallelism
instead.

## Usage

``` r
.MaximinMulti(k, firsts, d = NULL, points = NULL)
```

## Arguments

- k:

  Integer: target subsample size (`>= 1`).

- firsts:

  Integer vector of distinct first-selected indices, one per pass.

- d:

  Square pairwise distance matrix, or `NULL` when `points` is given.

- points:

  A `double` `N x dim` coordinate matrix, or `NULL`.

## Value

`.MaximinMulti()` returns a list with one `list(idx, tK)` per element of
`firsts`, in that order.
