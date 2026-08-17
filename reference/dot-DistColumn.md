# Normalise a distance-column oracle result to a masked length-`N` vector

The user's `colFn(i)` may either report the self-distance (a length-`N`
vector, position `i` ignored) or omit it (a length-`N - 1` vector of the
distances to the other elements, in index order). Either way this
returns a length-`N` numeric vector with position `i` set to `-Inf`, so
the downstream [`which.max()`](https://rdrr.io/r/base/which.min.html) /
[`pmin.int()`](https://rdrr.io/r/base/Extremes.html) never re-select
`i`. The mask invariant for the oracle path lives here, not in the
callers.

## Usage

``` r
.DistColumn(colFn, i, N)
```

## Arguments

- colFn:

  Column oracle; see
  [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md).

- i:

  Integer 1-based index whose distance column is requested.

- N:

  Integer element count.

## Value

`.DistColumn()` returns a numeric vector of length `N`, masked to `-Inf`
at position `i`.
