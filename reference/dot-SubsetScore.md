# Score a Gonzalez subset by its minimum (or mean) pairwise distance

Score a Gonzalez subset by its minimum (or mean) pairwise distance

## Usage

``` r
.SubsetScore(d, idx, objective = c("min_pairwise", "mean_pairwise"))
```

## Arguments

- d:

  Full pairwise distance matrix.

- idx:

  Integer indices of selected rows/cols.

- objective:

  `"min_pairwise"` (default; canonical Gonzalez T_k) or
  `"mean_pairwise"`.

## Value

Numeric scalar; `NA` if `length(idx) < 2`.
