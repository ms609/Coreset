# Normalise a distance-column oracle result for the DropAdd records

Reports the self-distance as 0, cf. `-Inf` in
[`.DistColumn()`](https://ms609.github.io/MaxMin/reference/dot-DistColumn.md)

## Usage

``` r
.DropAddColumn(colFn, i, N)
```

## Arguments

- colFn:

  Column oracle; see
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md).

- i:

  Integer 1-based index whose distance column is requested.

- N:

  Integer element count.

## Value

`.DropAddColumn()` returns a numeric vector of length `N` whose position
`i` is `0`, matching a distance matrix's diagonal.
