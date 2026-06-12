# Format and print MaxMin solver results

Terse summaries of the objects returned by the MaxMin solvers.

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

`print.MaxMin()` returns `x`, invisibly. It is called for its
side-effect of printing `format(x)` to the console. `format.MaxMin()`
returns a character string describing a `MaxMinSelection` (from
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) and
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md)); it
reports its size, the selected indices, the algorithm (and if applicable
strategy), and the achieved \\T_k\\. A `MaxMinExact` (from
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md))
object additionally states whether optimality was proven.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
print(FarFirst(5L, dist(pts)))
#> 5 elements (4 14 26 5 28) selected by farthest-first (best of 5 strategies, 3 tied: random_furthest1, random_furthest2, random_furthest3), each at distance >= 1.765
```
