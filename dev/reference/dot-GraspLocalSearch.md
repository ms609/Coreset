# Fast local search with the extended-improvement criterion.

Iterates 1-swap moves on critical elements (those participating in a
min-distance edge). A swap is accepted if it strictly increases d\*, or
if it preserves d\* while reducing the count of pairs at d\*.

## Usage

``` r
.GraspLocalSearch(d, sel)
```

## Arguments

- d:

  Square distance matrix.

- sel:

  Integer vector of size k.

## Value

`.GraspLocalSearch()` returns an improved integer vector of size k.
