# Exact discrete k-centre optimum on small instances

`ExactKCentre()` finds an optimal solution to the discrete *k*-centre
problem.

## Usage

``` r
ExactKCentre(k, d, maxSeconds = 60)

ExactKCenter(k, d, maxSeconds = 60)
```

## Arguments

- k:

  Integer specifying maximum number of centres to identify, from 1 to
  `nrow(d)`.

- d:

  `dist` object or a square symmetric numeric distance matrix.

- maxSeconds:

  Wall-clock budget in seconds for the whole search. If it expires
  before the optimum is proven, the smallest radius proven feasible so
  far is returned, with the attribute `proven = FALSE`.

## Value

`ExactKCentre()` returns an integer vector of length \\\le k\\ listing
the chosen centres in ascending order. It has class
`c("KCentreExact", "KCentreSelection")` and attributes:

- radius:

  The covering radius achieved; the proven optimum when `proven` is
  `TRUE`, otherwise an upper bound.

- proven:

  Logical: `TRUE` if optimality is certified.

- time_s:

  Wall-clock seconds elapsed.

- N, k:

  Instance size and centre budget.

It prints as a one-line summary and indexes a matrix or data frame
directly. The `"KCentreSelection"` superclass means
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
and any generic written for that class work here too.

The covering optimum may be attained by fewer than `k` centres (extra
centres never help once coverage is achieved); the result then has
length `< k` and the reported `radius` is still the proven *k*-centre
optimum. The problem is NP-hard, so this is an external ground-truth
reference for small instances, not a scalable method.

## Details

The optimum covering radius is the smallest threshold `r`, over the
achieved distinct distances, for which `k` centres can cover every point
within `r`.

Each probe solves a minimum-cardinality *set-cover* integer program with
the `highs` MILP backend, the covering constraints held as a sparse
matrix – the covering dual of
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)'s
node-packing program. The search is warm-started from the
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
(CDSh) radius, a proven feasible upper bound that caps the binary
search, then bisects downward to the smallest feasible radius.

## Progress bar

In interactive sessions, a progress indicator is shown. To toggle, set
`options("MaxMin.progress" = FALSE)` (or `TRUE`).

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) for
the fast near-optimal heuristic;
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
for the dual MMDP optimum.

## Examples

``` r
# \donttest{
if (requireNamespace("highs", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {
  set.seed(1)
  pts <- matrix(rnorm(40), ncol = 2)
  d <- dist(pts)
  ExactKCentre(3L, d)
}
#> 3 centres (3 11 15) by exact MILP (highs), proven optimal, covering radius = 1.385
# }
```
