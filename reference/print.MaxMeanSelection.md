# Format and print Max-Mean solver results

Terse one-line and detailed summaries of the objects returned by
[`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md).

## Usage

``` r
# S3 method for class 'MaxMeanSelection'
format(x, ...)

# S3 method for class 'MaxMeanSelection'
print(x, ...)

# S3 method for class 'MaxMeanSelection'
summary(object, ...)
```

## Arguments

- x:

  A `MaxMeanSelection` object returned by
  [`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md).

- ...:

  Ignored; present for S3 compatibility.

- object:

  A `MaxMeanSelection` object returned by
  [`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md).

## Value

`print.MaxMeanSelection()` returns `x`, invisibly.
`format.MaxMeanSelection()` returns a character string reporting the
selection size, the selected indices, and the achieved max-mean
objective \\f(S)\\.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`print.MaxEntropy`](https://ms609.github.io/MaxMin/reference/print.MaxEntropy.md),
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md),
[`print.MaxSum`](https://ms609.github.io/MaxMin/reference/print.MaxSum.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
print(MaxMean(dist(pts), maxSeconds = 1))
#> 30 elements (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ... (+10 more)) selected by MaxMean RLTS, f = 22.17
```
