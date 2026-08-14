# Exact Maximum Diversity Problem (max-sum) solution

`ExactMaxSum()` finds the optimal solution to the Max-Sum Diversity
Problem (the "maximum diversity problem"): select the `k`-subset of
points maximising the **total** pairwise distance it contains. It is the
max-sum counterpart of
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
(which maximises the *minimum* pairwise distance), solved by per-node
integer-program linearisation (Kuo et al. 1993) . As the problem is
NP-hard it is feasible only for small sets.

## Usage

``` r
ExactMaxSum(k, d, maxSeconds = 60, warmStart = NULL)
```

## Arguments

- k:

  Integer: target subset size, between 2 and `nrow(d)`.

- d:

  `dist` object or a square symmetric numeric distance matrix.

- maxSeconds:

  Numeric: search terminates after this many seconds have elapsed,
  returning largest threshold proven feasible.

- warmStart:

  Integer vector giving indices of a candidate subset to add to the
  heuristic warm-start pool, e.g. a selection computed by another
  solver.

## Value

`ExactMaxSum()` returns an integer vector of length `k` (sorted
ascending) with class `"MaxSumSelection"`, carrying attributes:

- score:

  Achieved total pairwise distance within the selection. When `proven`
  is `TRUE` this is the optimum; otherwise a lower bound.

- proven:

  Logical: `TRUE` if optimality was certified within the budget.

- time_s, N, k:

  Wall-clock seconds, instance size, target size.

## Details

The optimum is floored by a multi-start 1-swap local search, which warms
the lower bound and is returned when the MILP cannot prove optimality
within `maxSeconds` – so the result is always at least a strong
heuristic incumbent.

## References

Kuo C, Glover F, Dhir KS (1993). “Analyzing and modeling the maximum
diversity problem by zero-one programming.” *Decision Sciences*,
**24**(6), 1171–1185.
[doi:10.1111/j.1540-5915.1993.tb00509.x](https://doi.org/10.1111/j.1540-5915.1993.tb00509.x)
.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(20), ncol = 2)
ExactMaxSum(3L, dist(pts))
#> 3 elements (1 3 4) by exact MILP, proven optimal, total distance = 9.388
```
