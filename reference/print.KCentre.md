# Format and print k-centre solver results

Terse summaries of the objects returned by
[`KCentre()`](https://ms609.github.io/Coreset/reference/KCentre.md)
(`"KCentreSelection"`) and
[`ExactKCentre()`](https://ms609.github.io/Coreset/reference/ExactKCentre.md)
(`"KCentreExact"`)

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

`print.KCentre()` returns `x`, invisibly. It is called for its
side-effect of printing `format(x)` to the console. `format.KCentre()`
returns a character string describing a `KCentreSelection`; it reports
the centre count, the chosen indices, the method, and the achieved
covering radius (with proof status for the exact solver).

## See also

Other reporting functions:
[`print.Coreset`](https://ms609.github.io/Coreset/reference/print.Coreset.md),
[`print.MaxEntropy`](https://ms609.github.io/Coreset/reference/print.MaxEntropy.md),
[`print.MaxMeanSelection()`](https://ms609.github.io/Coreset/reference/print.MaxMeanSelection.md),
[`print.MaxSum`](https://ms609.github.io/Coreset/reference/print.MaxSum.md),
[`summary.Coreset`](https://ms609.github.io/Coreset/reference/summary.Coreset.md)

## Examples

``` r
set.seed(1)
KCentre(4L, dist(matrix(rnorm(60), ncol = 2)))
#> 4 centres (10 12 14 28) by CDSh, covering radius <= 1.242
```
