# Compose the one-line selection summary shared by both print methods

Compose the one-line selection summary shared by both print methods

## Usage

``` r
.MaxMinSummaryLine(n, idx, by, tk, maxShow = 20L)
```

## Arguments

- n:

  Integer count of selected elements.

- idx:

  Integer selected indices in stored order.

- by:

  Character phrase naming the algorithm (the "selected by ..." part).

- tk:

  Numeric achieved \\T_k\\ (the minimum pairwise distance), or `NA` when
  there is no pairwise distance to report (a single element, or all
  elements selected).

- maxShow:

  Integer index-list truncation threshold; see
  [`.FormatIndexList()`](https://ms609.github.io/MaxMin/reference/dot-FormatIndexList.md).

## Value

`.MaxMinSummaryLine()` returns a length-1 character string.
