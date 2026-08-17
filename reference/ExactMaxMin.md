# Exact Max-Min Diversity Problem solution

`ExactMaxMin()` finds the optimal solution to the Max-Min Diversity
Problem (discrete *p*-dispersion) by iterated node-packing (Sayyady and
Fathi 2016) (which may be slow or intractable on large sets).

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

  The minimum pairwise distance within the selection. When `proven` is
  `TRUE` this is the optimum; otherwise a lower bound.

- proven:

  Logical: `TRUE` if the search certified optimality within the budget,
  `FALSE` if it returned an unproven incumbent.

- time_s:

  Wall-clock seconds elapsed.

- N, k:

  Instance size and target subset size.

Prints as a terse summary via
[`print.MaxMinSelection()`](https://ms609.github.io/Coreset/reference/print.Coreset.md).

## Details

The search is warm-started from a heuristic lower bound (the best of
several [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)
restarts and a
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md)
pass), then gallops upward from that bound to the first infeasible
threshold and bisects the resulting bracket. When a heuristic already
attains the optimum, a single infeasibility proof certifies it. Each
feasibility probe is first reduced to its \\(k-1)\\-core and greedily
coloured, then searched exhaustively for a witness under a colouring
bound. The search runs on one core: its branches parallelise, but
measurably only for infeasibility proofs, and threads would make the
reported subset thread-dependent. The indices returned may vary between
releases where several subsets attain the optimum; the `score` does not.

## Progress bar

In interactive sessions, a progress indicator is shown. To toggle, set
`options("Coreset.progress" = FALSE)` (or `TRUE`).

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
#> 3 elements (4 5 6) selected by exact solver, proven optimal, each at distance >= 2.035
```
