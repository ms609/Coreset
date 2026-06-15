# GRASP with Path Relinking for the Max-Min Diversity Problem

`Grasp()` solves the Max-Min Diversity Problem (discrete *p*-dispersion)
with the static variant of the GRASP / path-relinking metaheuristic
(Resende et al. 2010, fig. 4) . This is the most expensive heuristic in
this package, and attains correspondingly high-quality selections.

## Usage

``` r
Grasp(k, d, plateau = 100L, eliteSize = 10L, alpha = 0.8, maxSeconds = Inf)
```

## Arguments

- k:

  Integer subset size, `2 <= k <= nrow(d)`.

- d:

  Either a `dist` object or a square symmetric numeric matrix.

- plateau:

  Integer; stop after this many consecutive GRASP iterations without an
  improvement to the best elite objective. The primary, deterministic
  stopping criterion. Default 100.

- eliteSize:

  Size of the elite set \|ES\|. Default 10.

- alpha:

  RCL threshold; `alpha = 1` is pure greedy, `alpha = 0` uniform random.
  Default 0.8.

- maxSeconds:

  Optional wall-clock ceiling in seconds. Default `Inf` (no ceiling,
  fully reproducible). A finite value caps runtime but makes the result
  machine-dependent.

## Value

`Grasp()` returns an integer vector of length `k` (1-based) **sorted
ascending** (unlike
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
which returns farthest-first order) with attributes:

- score:

  Achieved MaxMin objective \\T_k\\.

- time_s:

  Wall-clock seconds spent.

- iters:

  Number of GRASP refinement iterations executed.

- pr_calls:

  Number of path-relinking pair-applications run.

The vector has class `"MaxMinSelection"` and prints as a one-line
summary (see
[`print.MaxMinSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md));
it is otherwise an ordinary integer vector.

## Details

The GRASP with path-relinking algorithm conducts a randomised-greedy
construction with extended-improvement local search builds; it
identifies an elite set, then conducts a single pass of path relinking
over all elite pairs (Resende et al. 2010) .

The refinement loop stops after `plateau` consecutive GRASP iterations
fail to improve the best elite objective, or once `maxSeconds` have
elapsed.

This method will fail if the complete \\N \times N\\ distance matrix, is
too large to fit into memory.

## Progress bar

In interactive sessions, a bar tracks how close the search is to its
`plateau` stopping criterion, snapping back each time a better solution
is found. To toggle, set `options("MaxMin.progress" = FALSE)` (or
`TRUE`).

## References

Resende MGC, Martí R, Gallego M, Duarte A (2010). “GRASP and path
relinking for the max-min diversity problem.” *Computers & Operations
Research*, **37**(3), 498–508.
[doi:10.1016/j.cor.2008.05.011](https://doi.org/10.1016/j.cor.2008.05.011)
.

## See also

[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) for
scalable refinement;
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
for the proven optimum on small instances.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
Grasp(5L, dist(pts), plateau = 20L, eliteSize = 4L)
#> 5 elements (3 4 5 24 25) selected by GRASP with path-relinking, each at distance >= 1.778
```
