# Format and print k-centre solver results

One-line summaries of the objects returned by
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
(`"KCentreSelection"`) and
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
(`"KCentreExact"`): the centre count, the chosen indices, the method,
and the achieved covering radius (with proof status for the exact
solver). Both objects are otherwise unchanged.

## Usage

``` r
# S3 method for class 'KCentreSelection'
format(x, ...)

# S3 method for class 'KCentreSelection'
print(x, ...)

# S3 method for class 'KCentreExact'
format(x, ...)

# S3 method for class 'KCentreExact'
print(x, ...)
```

## Arguments

- x:

  A `"KCentreSelection"` or `"KCentreExact"` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

`print.KCentre()` returns `x`, invisibly (`print`); a length-1 character
string (`format`).

## See also

Other reporting functions:
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
KCentre(4L, dist(matrix(rnorm(60), ncol = 2)))
#> 4 centres (10 12 14 28) by CDSh, covering radius <= 1.242
```
