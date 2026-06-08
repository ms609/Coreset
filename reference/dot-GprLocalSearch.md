# Fast local search with the extended-improvement criterion.

Iterates 1-swap moves on critical elements (those participating in a
min-distance edge). A swap is accepted if it strictly increases d\*, or
if it preserves d\* while reducing the count of pairs at d\*.

## Usage

``` r
.GprLocalSearch(d, sel)
```

## Arguments

- d:

  Square distance matrix.

- sel:

  Integer vector of size m.

## Value

Improved integer vector of size m.
