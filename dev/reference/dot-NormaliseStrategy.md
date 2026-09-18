# Read start indices out of a character `strategy`

[`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md)
takes a start index wherever it takes a strategy name, so an index can
share a vector with names: R coerces `c(17, "random_furthest")` to
`c("17", "random_furthest")`. Each element that reads as a number is
rewritten in canonical integer form (`"1e+05"` becomes `"100000"`), and
`"first"` is a synonym for `"1"`. Names are left for the caller to
validate.

## Usage

``` r
.NormaliseStrategy(strategy)
```

## Arguments

- strategy:

  A character (or multi-element numeric) `strategy`.

## Value

`.NormaliseStrategy()` returns `strategy` as a character vector with
every start index written as a positive integer string, or stops.
