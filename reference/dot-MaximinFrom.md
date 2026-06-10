# Gonzalez maximin from a single starting index

Internal helper: greedy furthest-point selection starting from a
specified index.

## Usage

``` r
.MaximinFrom(d, m, first)
```

## Arguments

- d:

  Square pairwise distance matrix.

- m:

  Integer: target subsample size (`>= 1`).

- first:

  Integer: index of the first selected point.

## Value

Integer vector of length `m` of selected row/col indices.
