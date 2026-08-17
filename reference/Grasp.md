# GRASP with Path Relinking for the Max-Min Diversity Problem

`Grasp()` solves the Max-Min Diversity Problem (discrete *p*-dispersion)
with the static variant of the GRASP / path-relinking metaheuristic
(Resende et al. 2010, fig. 4) . This is the most expensive heuristic in
this package, and attains correspondingly high-quality selections.

## Usage

``` r
Grasp(
  k,
  d,
  plateau = 100L,
  eliteSize = 10L,
  alpha = 0.8,
  maxSeconds = Inf,
  maxCandidates = 2000L
)
```

## Arguments

- k:

  Integer subset size, `2 <= k <= nrow(d)`.

- d:

  Either a `dist` object or a square symmetric numeric matrix.

- plateau:

  Integer; stop after this many consecutive GRASP iterations without an
  improvement to the best elite objective. The primary, deterministic
  stopping criterion.

- eliteSize:

  Size of the elite set \|ES\|.

- alpha:

  Construction greediness in `[0, 1]`. Each step draws the next point at
  random from a shortlist of the strongest candidates – those whose gain
  lies within a fraction `alpha` of the best-to-worst spread.
  `alpha = 1` is pure greedy (best only); `alpha = 0` is uniform random
  among candidates.

- maxSeconds:

  Numeric specifying wall-clock ceiling, in seconds.

- maxCandidates:

  Integer: a composable-coreset tractability cap. When the number of
  candidate points `N` exceeds `maxCandidates`,
  [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
  thins the candidates to a `maxCandidates`-point coreset (with the
  deterministic, RNG-free `"peripheral"` seed, so no random stream is
  perturbed), the solver runs on the coreset, and the chosen indices are
  mapped back to the original numbering. This lets the solver produce a
  solution at scales where it would otherwise be intractable.
  `maxCandidates = 0` (or `Inf`) disables thinning and runs on the full
  problem; a cap at or above `N` is a no-op. A cap below `k` is an
  error. The default is `2000L`, conservative because `Grasp()` is
  matrix-only, so Thinning is **on by default**: an input larger than
  the cap is thinned (and a warning is emitted) unless
  `maxCandidates = 0` is passed.

## Value

`Grasp()` returns an integer vector of length `k` specifying the indices
of the selected points, with attributes:

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
[`print.MaxMinSelection()`](https://ms609.github.io/Coreset/reference/print.Coreset.md)).

## Details

The GRASP with path-relinking algorithm conducts a randomised-greedy
construction with extended-improvement local search builds; it
identifies an elite set, then conducts a single pass of path relinking
over all elite pairs (Resende et al. 2010) .

The refinement loop stops after `plateau` consecutive GRASP iterations
fail to improve the best elite objective, or once `maxSeconds` have
elapsed.

This method will fail if the complete \\N \times N\\ distance matrix is
too large to fit into memory.

## Parallelism

To parallelize computation when OpenMP is available, set the
`"mc.cores"` option:

    options(mc.cores = 2L)                       # use a fixed number of cores
    options(mc.cores = parallel::detectCores())  # or all available cores

## Progress bar

In interactive sessions, a bar tracks how close the search is to its
`plateau` stopping criterion, snapping back each time a better solution
is found. To toggle, set `options("Coreset.progress" = FALSE)` (or
`TRUE`).

## References

Resende MGC, Martí R, Gallego M, Duarte A (2010). “GRASP and path
relinking for the max-min diversity problem.” *Computers & Operations
Research*, **37**(3), 498–508.
[doi:10.1016/j.cor.2008.05.011](https://doi.org/10.1016/j.cor.2008.05.011)
.

## See also

[`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md) for
scalable refinement;
[`ExactMaxMin()`](https://ms609.github.io/Coreset/reference/ExactMaxMin.md)
for the proven optimum on small instances.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
Grasp(5L, dist(pts), plateau = 20L, eliteSize = 4L)
#> 5 elements (3 4 5 24 25) selected by GRASP with path-relinking, each at distance >= 1.778

# Composable coreset: thin to 20 candidates with farthest-first, then run
# GRASP on the coreset. Returned indices are original-space row indices.
suppressWarnings(Grasp(5L, dist(pts), plateau = 20L, maxCandidates = 20L))
#> 5 elements (3 4 5 24 25) selected by GRASP with path-relinking, each at distance >= 1.778
```
