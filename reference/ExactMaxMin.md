# Exact Max-Min Diversity Problem optimum on small instances

Solves the Max-Min Diversity Problem (discrete p-dispersion) to proven
optimality by iterated node-packing Sayyady and Fathi (2016) : the
optimum is the largest threshold `lambda`, over the achieved distinct
pairwise distances, for which the threshold graph `G(lambda)` (edges
join pairs closer than `lambda`) contains an independent set of size at
least `k`. Each probe solves a maximum-independent-set integer program
with the `highs` MILP backend, the packing constraints held as a sparse
matrix.

## Usage

``` r
ExactMaxMin(k, d, solver = NULL, maxSeconds = 60, warmStart = NULL)
```

## Arguments

- k:

  Integer target subset size, `2 <= k <= nrow(d)`.

- d:

  A `dist` object or a square symmetric numeric distance matrix.

- solver:

  Solver to use. Currently only `"highs"` is implemented; `NULL` selects
  it. Other values raise an error.

- maxSeconds:

  Wall-clock budget in seconds for the whole search (shared across all
  internal IP solves). If the budget expires before the optimum is
  proven, the largest threshold proven feasible so far is returned with
  `proven = FALSE`.

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
one-line summary (size, indices, solver, proof status and achieved
`T_k`; see
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
attains the optimum – common at the small `k` for which an exact
reference is wanted – a single infeasibility solve certifies it. The
warm start only sets the starting lower bound: the returned optimum is
proven regardless of heuristic quality (a loose seed merely costs extra
solves). The problem is NP-hard, so this remains an external
ground-truth reference for small instances, not a scalable method.

The proven `objective` is exact and does not depend on the RNG. Only the
returned `indices` can vary when several subsets attain the optimum: the
warm start draws on the session RNG via
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) and, like
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md), advances
it. Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before
`ExactMaxMin()` for a reproducible selection.

## Progress bar

Shows a progress indicator controlled by
`getOption("MaxMin.progress", interactive())` — `TRUE` by default in
interactive sessions, `FALSE` otherwise.

## References

Sayyady F, Fathi Y (2016). “An integer programming approach for solving
the p-dispersion problem.” *European Journal of Operational Research*,
**253**(1), 216–225.
[doi:10.1016/j.ejor.2016.02.026](https://doi.org/10.1016/j.ejor.2016.02.026)
.
