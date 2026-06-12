# Gonzalez maximin from a distance-column oracle (worker)

Mirrors `MaximinFrom_cpp()`, substituting an on-demand `colFn(i)` call
for the matrix-column read `d[, i]`.
[`which.max()`](https://rdrr.io/r/base/which.min.html) uses R's
first-maximum (strict `>`) rule, matching the kernel's tie-breaking, so
the selection is identical to the matrix path on symmetric input.

## Usage

``` r
.MaximinFromColumn(colFn, N, k, first, progress = FALSE)
```

## Arguments

- colFn:

  Column oracle; see
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md).

- N:

  Integer element count.

- k:

  Integer subset size (`>= 2`).

- first:

  Integer seed index.

## Value

Integer vector of selected indices.
