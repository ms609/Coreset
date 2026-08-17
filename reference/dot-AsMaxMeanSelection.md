# Stamp the `MaxMeanSelection` class onto a [`MaxMean()`](https://ms609.github.io/Coreset/reference/MaxMean.md) result

Parallel to
[`.AsMaxMinSelection()`](https://ms609.github.io/Coreset/reference/dot-AsMaxMinSelection.md)
for the fixed-cardinality solvers. An empty selection is returned
unchanged.

## Usage

``` r
.AsMaxMeanSelection(x)
```

## Arguments

- x:

  Integer index vector carrying `score`, `size`, `time_s`, `iters`.

## Value

`.AsMaxMeanSelection()` returns `x` with class `"MaxMeanSelection"`, or
`x` unchanged if it is empty.
