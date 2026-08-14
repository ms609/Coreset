# Format and print Max-Sum (maximum diversity) solver results

Terse summary of the object returned by
[`ExactMaxSum()`](https://ms609.github.io/MaxMin/reference/ExactMaxSum.md)
(`"MaxSumSelection"`), reporting the achieved **total** pairwise
distance (the max-sum objective) rather than the minimum distance of
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md).

## Usage

``` r
# S3 method for class 'MaxSumSelection'
format(x, ...)

# S3 method for class 'MaxSumSelection'
print(x, ...)
```

## Arguments

- x:

  A `"MaxSumSelection"` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

`format.MaxSumSelection()` returns a one-line character summary;
`print.MaxSumSelection()` returns `x` invisibly, called for its
side-effect.

## See also

Other reporting functions:
[`print.KCentre`](https://ms609.github.io/MaxMin/reference/print.KCentre.md),
[`print.MaxEntropy`](https://ms609.github.io/MaxMin/reference/print.MaxEntropy.md),
[`print.MaxMeanSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMeanSelection.md),
[`print.MaxMin`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md),
[`summary.MaxMin`](https://ms609.github.io/MaxMin/reference/summary.MaxMin.md)

## Examples

``` r
set.seed(1)
ExactMaxSum(3L, dist(matrix(rnorm(20), ncol = 2)))
#> 3 elements (1 3 4) by exact MILP, proven optimal, total distance = 9.388
```
