# Format and print MaxMin solver results

One-line, human-readable summaries of the objects returned by the MaxMin
solvers. A `MaxMinSelection` (from
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) and
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md)) reports
its size, the selected indices, the algorithm (and, for a
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
ensemble, the winning strategy), and the achieved \\T_k\\; a
`MaxMinExact` (from
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md))
additionally states whether optimality was proven. `MaxMinExact` extends
`MaxMinSelection` (`class = c("MaxMinExact", "MaxMinSelection")`), so
`inherits(x, "MaxMinSelection")` is `TRUE` for any solver result. Both
objects are otherwise unchanged – a `MaxMinSelection` still indexes like
the bare integer vector it wraps – so these methods only affect display.

## Usage

``` r
# S3 method for class 'MaxMinSelection'
format(x, ...)

# S3 method for class 'MaxMinSelection'
print(x, ...)

# S3 method for class 'MaxMinExact'
format(x, ...)

# S3 method for class 'MaxMinExact'
print(x, ...)
```

## Arguments

- x:

  A `MaxMinSelection` or `MaxMinExact` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

`x`, invisibly (`print`); a length-1 character string (`format`).

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
FarFirst(dist(pts), 5L)
#> 5 elements (14 4 26 5 28) selected by Gonzalez farthest-first (best of 3 strategies, winner random_furthest2), each at distance >= 1.765
```
