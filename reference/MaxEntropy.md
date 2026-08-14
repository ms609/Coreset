# Maximum-entropy (maxdet) subset selection

`MaxEntropy()` selects the `k`-subset of points that maximises the
log-determinant of its kernel block, \\\log\det K_S\\ – the spanned
volume of the selection, the maximum-entropy sampling criterion (Shewry
and Wynn 1987) and the maximum-a-posteriori mode of a determinantal
point process (Kulesza and Taskar 2012) . A redundant point lies in the
span of those already chosen, adds zero volume, and is never taken, so
the selection is exactly density-blind.

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

  Integer target selection size, \\1 \le k \le n\\.

- d:

  A `dist` object or square numeric distance matrix over the `n` points.

- sigma:

  Optional kernel bandwidth; defaults to the median positive distance.

- repair:

  Positive Semi-Definite repair method for the kernel: `"clip"`
  (nearest), `"shift"` (diagonal loading) or `"truncate"` (low-rank
  embedding).

- exact:

  Logical, or `NA` (the default) to choose automatically: use the exact
  enumeration when `choose(n, k) <= maxCombos`, otherwise the greedy.
  `TRUE` forces enumeration (error if it exceeds `maxCombos`); `FALSE`
  forces the greedy.

- maxCombos:

  Numeric ceiling on `choose(n, k)` for exact enumeration.

## Value

`MaxEntropy()` returns an integer vector of length `k` (sorted
ascending) with class `"MaxEntropySelection"`, carrying attributes:

- logDet, score:

  The retained \\\log\det K_S\\ of the selection; `-Inf` for a
  degenerate selection (one forced to repeat near-identical points
  because `k` exceeds the number of distinct points, which also warns).

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
built from the supplied distances (`sigma` defaulting to the median
positive distance) and repaired to a positive-semidefinite matrix,
because a general distance is not of negative type and \\\log\det\\
requires it. The exact argmax is NP-hard (Kulesza and Taskar 2012) , so
the selection is built greedily by pivoted Cholesky – at each step
adding the point of largest residual conditional variance – with exact
enumeration substituted where \\\binom{n}{k}\\ does not exceed
`maxCombos`. The greedy first pivot is tied on a unit-diagonal kernel
and is broken deterministically by the most peripheral point (least
total similarity), so no random seed is used.

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
