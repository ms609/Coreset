# Exact Max-Min Diversity (MMDP) optimum on small instances

Solves the Max-Min Diversity Problem (discrete p-dispersion) to proven
optimality by iterated node-packing Sayyady and Fathi (2016) : the
optimum is the largest threshold `lambda`, over the achieved distinct
pairwise distances, for which the threshold graph `G(lambda)` (edges
join pairs closer than `lambda`) contains an independent set of size at
least `m`. A binary search over the sorted distances resolves that
threshold; each probe solves a maximum-independent-set integer program
with the `highs` MILP backend. The problem is NP-hard, so this is
intended only as an external ground-truth reference on small instances,
not a scalable method.

## Usage

``` r
ExactMaxMin(
  d,
  m,
  solver = NULL,
  time_budget_s = 60,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- d:

  A `dist` object or a square symmetric numeric distance matrix.

- m:

  Integer target subset size, `2 <= m <= nrow(d)`.

- solver:

  Solver to use. Currently only `"highs"` is implemented; `NULL` selects
  it. Other values raise an error.

- time_budget_s:

  Wall-clock budget in seconds for the whole search (shared across all
  internal IP solves). If the budget expires before the optimum is
  proven, the largest threshold proven feasible so far is returned with
  `proven = FALSE`.

- progress:

  Logical; show a progress bar during the binary search. Default: `TRUE`
  in interactive sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`).

## Value

A list with fields

- indices:

  Integer vector of length `m`, sorted ascending: the selected points.

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

- n, m:

  Instance size and target subset size.

## References

Sayyady F, Fathi Y (2016). “An integer programming approach for solving
the p-dispersion problem.” *European Journal of Operational Research*,
**253**(1), 216–225.
[doi:10.1016/j.ejor.2016.02.026](https://doi.org/10.1016/j.ejor.2016.02.026)
.
