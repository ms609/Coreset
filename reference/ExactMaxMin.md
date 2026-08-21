# Exact Max-Min Diversity Problem solution

`ExactMaxMin()` finds the optimal solution to the Max-Min Diversity
Problem (discrete *p*-dispersion) by iterated node-packing (Sayyady and
Fathi 2016) (which may be slow or intractable on large sets).

## Usage

``` r
ExactMaxMin(
  k,
  d,
  maxSeconds = 60,
  warmStart = NULL,
  nStart = 1L,
  graspPlateau = 50L,
  dropPlateau = 512L
)
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

  Optional integer vector giving indices of a candidate subset to add to
  the heuristic warm-start pool.

- nStart:

  Integer: how many
  [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)
  restarts enter the warm-start pool.

- graspPlateau, dropPlateau:

  Integer: the stopping plateaus given to the pool's
  [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)
  restarts and its
  [`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md)
  pass. Deeper searches cost more, but raise the lower bound the exact
  search starts from.

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
`nStart` [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)
restarts and a
[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md)
pass), then gallops upward from that bound to the first infeasible
threshold and bisects the resulting bracket.

To parallelize computation when OpenMP is available, set the
`"mc.cores"` option:


    options(mc.cores = 2L)                       # use a fixed number of cores
    options(mc.cores = parallel::detectCores())  # or all available cores

Parallelization returns identical results under a given seed.

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
