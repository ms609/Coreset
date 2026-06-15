# Exact Max-Min Diversity Problem solution

`ExactMaxMin()` finds the optimal solution to the Max-Min Diversity
Problem (discrete *p*-dispersion) by iterated node-packing (Sayyady and
Fathi 2016) . As this problem is NP-hard, it is feasible only for small
sets.

## Usage

``` r
ExactMaxMin(k, d, maxSeconds = 60, warmStart = NULL)
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

`ExactMaxMin()` returns an integer vector of length `k` (sorted
ascending) with class `"MaxMinSelection"`, carrying attributes:

- score:

  Achieved `T_k` – the minimum pairwise distance within the selection.
  When `proven` is `TRUE` this is the true optimum; otherwise a lower
  bound.

- proven:

  Logical: `TRUE` if the search certified optimality within the budget,
  `FALSE` if it returned an unproven incumbent.

- time_s:

  Wall-clock seconds elapsed.

- N, k:

  Instance size and target subset size.

Prints as a one-line summary via
[`print.MaxMinSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md).
`inherits(result, "MaxMinSelection")` is `TRUE`.

## Details

The search is warm-started from a heuristic lower bound (the best of
several [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md)
restarts and a
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md)
pass), then gallops upward from that bound to the first infeasible
threshold and bisects the resulting bracket. When a heuristic already
attains the optimum, a single infeasibility solve certifies it.

The proven `objective` is exact and does not depend on the RNG. Only the
returned `indices` can vary when several subsets attain the optimum: the
warm start draws on the session RNG via
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) and, like
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md), advances
it. Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before
`ExactMaxMin()` for a reproducible selection.

## Progress bar

In interactive sessions, a progress indicator is shown. To toggle, set
`options("MaxMin.progress" = FALSE)` (or `TRUE`).

## References

Sayyady F, Fathi Y (2016). “An integer programming approach for solving
the p-dispersion problem.” *European Journal of Operational Research*,
**253**(1), 216–225.
[doi:10.1016/j.ejor.2016.02.026](https://doi.org/10.1016/j.ejor.2016.02.026)
.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(18), ncol = 2)
ExactMaxMin(3L, dist(pts))
#> 3 elements (3 4 5) selected by exact solver, proven optimal, each at distance >= 2.035
```
