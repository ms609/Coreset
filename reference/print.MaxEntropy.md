# Format and print maximum-entropy (maxdet) solver results

Terse summary of the object returned by
[`MaxEntropy()`](https://ms609.github.io/MaxMin/reference/MaxEntropy.md)
(`"MaxEntropySelection"`), reporting the retained log-determinant
(\\\log\det K_S\\, the maxdet objective) and the magnitude of the
positive-semidefinite repair.

## Usage

``` r
# S3 method for class 'MaxEntropySelection'
format(x, ...)

# S3 method for class 'MaxEntropySelection'
print(x, ...)
```

## Arguments

- x:

  A `"MaxEntropySelection"` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

`format.MaxEntropySelection()` returns a one-line character summary;
`print.MaxEntropySelection()` returns `x` invisibly, called for its
side-effect.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`print.MaxMeanSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMeanSelection.md),
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md),
[`print.MaxSum`](https://ms609.github.io/MaxMin/reference/print.MaxSum.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
MaxEntropy(4L, dist(matrix(rnorm(40), ncol = 2)))
#> 4 elements (4 11 14 16) by max-entropy, exact enumeration, log det = -0.1964 (repair removed 0 of mass)
```
