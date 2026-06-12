# Greedy path relinking from x toward y.

Greedy path relinking from x toward y.

## Usage

``` r
.GraspPathRelink(d, x, y)
```

## Arguments

- d:

  Square distance matrix.

- x, y:

  Integer selections of equal length.

## Value

`.GraspPathRelink()` returns
`list(best = best selection on path, intermediates = number of intermediate states visited including endpoints)`.
