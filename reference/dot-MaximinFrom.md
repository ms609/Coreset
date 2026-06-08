# Gonzalez maximin from a single starting index

Internal helper: greedy furthest-point selection starting from a
specified index.

## Usage

``` r
.MaximinFrom(d, n, first)
```

## Arguments

- d:

  Square pairwise distance matrix.

- n:

  Integer: target subsample size (`>= 1`).

- first:

  Integer: index of the first selected point.

## Value

Integer vector of length `n` of selected row/col indices.
