# Maximum-entropy (maxdet) subset selection

`MaxEntropy()` selects the `k` points that maximise the log-determinant
of their kernel block, \\\log\det K_S\\. This corresponds to the volume
spanned by the selection, the maximum-entropy sampling criterion (Shewry
and Wynn 1987) and the maximum-a-posteriori mode of a determinantal
point process (Kulesza and Taskar 2012) .

## Usage

``` r
MaxEntropy(
  k,
  d,
  sigma = NULL,
  repair = c("clip", "shift", "truncate"),
  exact = NA,
  maxCombos = 3e+05
)
```

## Arguments

- k:

  Integer specifying target selection size, \\1 \le k \le n\\.

- d:

  `dist` object or square numeric distance matrix over the `n` points.

- sigma:

  Optional numperic specifying kernel bandwidth; defaults to the median
  positive distance.

- repair:

  Character specifying positive Semi-Definite repair method for the
  kernel: `"clip"` (nearest), `"shift"` (diagonal loading) or
  `"truncate"` (low-rank embedding).

- exact:

  Logical: `TRUE` uses explicit enumeration, failing with an error if
  `maxCombos` is exceeded; `FALSE` uses the greedy approximation. `NA`
  uses exact enumeration when `choose(n, k) <= maxCombos`, greedy
  otherwise.

- maxCombos:

  Numeric specifying ceiling on `choose(n, k)` for exact enumeration.

## Value

`MaxEntropy()` returns an integer vector of length `k` (sorted
ascending) with class `"MaxEntropySelection"`, carrying attributes:

- logDet, score:

  The retained \\\log\det K_S\\ of the selection. `-Inf` is returned for
  a degenerate selection where `k` exceeds the number of distinct
  points.

- negMass:

  Fraction of spectral mass removed by the Positive Semi-Definite
  repair.

- sigma, repair, exact:

  The bandwidth, repair, and whether the optimum was certified by
  enumeration.

- seed, N, k:

  The peripheral seed index, instance size, target size.

## Details

A radial-basis kernel \\K\_{ij} = \exp(-d\_{ij}^2 / 2\sigma^2)\\ is
built from the supplied distances and repaired to a
positive-semidefinite matrix. The exact argmax is NP-hard (Kulesza and
Taskar 2012) . A greedy approximation is built by pivoted Cholesky,
adding at each step the point of largest residual conditional variance.
Ties are broken by selecting the more peripheral point.

## References

Kulesza A, Taskar B (2012). “Determinantal point processes for machine
learning.” *Foundations and Trends in Machine Learning*, **5**(2–3),
123–286. [doi:10.1561/2200000044](https://doi.org/10.1561/2200000044)
.  
  
Shewry MC, Wynn HP (1987). “Maximum entropy sampling.” *Journal of
Applied Statistics*, **14**(2), 165–170.
[doi:10.1080/02664768700000020](https://doi.org/10.1080/02664768700000020)
.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(40), ncol = 2)
MaxEntropy(4L, dist(pts))
#> 4 elements (4 11 14 16) by max-entropy, exact enumeration, log det = -0.1964 (repair removed 0 of mass)
```
