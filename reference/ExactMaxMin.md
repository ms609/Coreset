# Exact Max-Min Diversity Problem solution

`ExactMaxMin()` finds the optimal solution to the Max-Min Diversity
Problem (discrete *p*-dispersion) by iterated node-packing Sayyady and
Fathi (2016) . As this problem is NP-hard, it is feasible only for small
sets.

## Usage

``` r
ExactMaxMin(k, d, maxSeconds = 60, warmStart = NULL)
```

## Arguments

- k:

  Integer target subset size, `2 <= k <= nrow(d)`.

- d:

  A `dist` object or a square symmetric numeric distance matrix.

- maxSeconds:

  Numeric: search terminates after this many seconds have elapsed;
  returning largest threshold proven feasible so far.

- warmStart:

  Optional integer vector: a candidate `k`-subset (1-based indices into
  `d`) to add to the heuristic warm-start pool, e.g. a selection already
  computed by another solver. Ignored unless it is a valid `k`-subset.
  The internal heuristics run regardless; a good `warmStart` can only
  reduce the number of IP solves, never change the proven optimum.

## Value

`ExactMaxMin()` returns a list (unlike
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) and
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
which each return a bare integer vector carrying a `score` attribute),
since it reports both the optimum and a proof status. The fields are

- indices:

  Integer vector of length `k`, sorted ascending: the selected points.

- objective:

  The achieved `T_k` – the minimum pairwise distance within `indices`.
  When `proven` is `TRUE` this equals the threshold `lambda` and is the
  true optimum; otherwise it is a valid lower bound on the optimum.

- proven:

  Logical: `TRUE` if the search certified optimality within the budget,
  `FALSE` if it returned an unproven incumbent.

- time_s:

  Wall-clock seconds elapsed.

- solver:

  Name of the MILP backend used.

- n, k:

  Instance size and target subset size.

The list has class `c("MaxMinExact", "MaxMinSelection")` and prints as a
one-line summary (size, indices, proof status and achieved `T_k`; see
[print.MaxMinSelection](https://ms609.github.io/MaxMin/reference/print.MaxMin.md));
it is otherwise an ordinary list. The `"MaxMinSelection"` superclass
means `inherits(result, "MaxMinSelection")` is `TRUE`, and any generic
written against that class works here too.

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
