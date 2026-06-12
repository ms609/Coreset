# Exact discrete k-centre optimum on small instances

Solves the discrete (vertex) *k*-centre problem to proven optimality:
the optimum covering radius is the smallest threshold `r`, over the
achieved distinct distances, for which `k` centres can cover every point
within `r`. Each probe solves a minimum-cardinality *set-cover* integer
program with the `highs` MILP backend, the covering constraints held as
a sparse matrix – the covering dual of
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)'s
node-packing program. The search is warm-started from the
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
(CDSh) radius, a proven feasible upper bound that caps the binary
search, then bisects downward to the smallest feasible radius.

## Usage

``` r
ExactKCentre(
  k,
  d,
  solver = NULL,
  maxSeconds = 60,
  warmStart = NULL,
  progress = getOption("MaxMin.progress", interactive())
)

ExactKCenter(
  k,
  d,
  solver = NULL,
  maxSeconds = 60,
  warmStart = NULL,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- k:

  Integer centre budget, `1 <= k <= nrow(d)`.

- d:

  A `dist` object or a square symmetric numeric distance matrix.

- solver:

  Solver to use. Currently only `"highs"` is implemented; `NULL` selects
  it.

- maxSeconds:

  Wall-clock budget in seconds for the whole search (shared across the
  internal IP solves). If it expires before the optimum is proven, the
  smallest radius proven feasible so far is returned with
  `proven = FALSE`.

- warmStart:

  Currently unused; reserved for a caller-supplied feasible centre set.
  The internal CDSh warm start runs regardless.

- progress:

  Logical; show a progress indicator. Default: `TRUE` in interactive
  sessions (`getOption("MaxMin.progress", interactive())`).

## Value

`ExactKCentre()` returns a list of class `"KCentreExact"` with fields

- indices:

  Integer vector (ascending), length `<= k`: the centres.

- radius:

  The covering radius they achieve; the proven optimum when `proven` is
  `TRUE`, otherwise a valid upper bound.

- proven:

  Logical: `TRUE` if optimality was certified within budget.

- time_s:

  Wall-clock seconds elapsed.

- solver:

  Name of the MILP backend used.

- n, k:

  Instance size and centre budget.

- n_centres:

  `length(indices)`.

It prints as a one-line summary and is otherwise an ordinary list.

## Details

The covering optimum may be attained by fewer than `k` centres (extra
centres never help once coverage is achieved); `indices` then has length
`< k` and the reported `radius` is still the proven *k*-centre optimum.
The problem is NP-hard, so this is an external ground-truth reference
for small instances, not a scalable method.

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
