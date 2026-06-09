# GRASP with Path Relinking for the Max-Min Diversity Problem

Solves the Max-Min Diversity Problem (discrete p-dispersion) with the
GRASP / path-relinking metaheuristic of Resende et al. (2010) , static
variant (their Fig. 4): a randomised-greedy construction with
extended-improvement local search builds and maintains an elite set,
followed by a single pass of path relinking over all elite pairs. On the
application benchmark this attains the highest \\T_k\\ of the methods in
this package, at correspondingly higher cost.

## Usage

``` r
GraspPR(
  d,
  m,
  max_no_improve = 100L,
  max_iter = NULL,
  elite_size = 10L,
  alpha = 0.8,
  time_budget_s = Inf,
  seed = NULL
)
```

## Arguments

- d:

  Either a `dist` object or a square symmetric numeric matrix.

- m:

  Integer subset size, `2 <= m <= nrow(d)`.

- max_no_improve:

  Integer; stop after this many consecutive GRASP iterations without an
  improvement to the best elite objective. The primary, deterministic
  stopping criterion. Default 100.

- max_iter:

  Optional integer hard cap on GRASP refinement iterations (excluding
  the elite-set construction). `NULL` (default) leaves `max_no_improve`
  in sole control.

- elite_size:

  Size of the elite set \|ES\|. Default 10.

- alpha:

  RCL threshold; `alpha = 1` is pure greedy, `alpha = 0` uniform random.
  Default 0.8.

- time_budget_s:

  Optional wall-clock ceiling in seconds. Default `Inf` (no ceiling,
  fully reproducible). A finite value caps runtime but makes the result
  machine-dependent.

- seed:

  Optional integer; if supplied, `set.seed(seed)` is called at entry.
  `GraspPR` is genuinely stochastic (randomised construction and RCL
  sampling), so the seed governs the trajectory and the returned
  selection.

## Value

A list with elements

- indices:

  Integer vector of length `m`, 1-based, sorted ascending.

- objective:

  Achieved MaxMin objective \\T_k\\.

- time_s:

  Wall-clock seconds spent.

- iters:

  Number of GRASP refinement iterations executed.

- pr_calls:

  Number of path-relinking pair-applications run.

## Details

**Deterministic termination.** The refinement loop stops after
`max_no_improve` consecutive GRASP iterations that fail to improve the
best elite objective (rather than after a wall-clock budget). Given a
fixed `seed`, the entire run — construction RNG, iteration count, and
result — is therefore reproducible and machine-independent;
`set.seed(seed)` controls the same random stream the compiled kernel
consumes via R's RNG. An optional `time_budget_s` ceiling is available
as a safety cap, but using a finite value reintroduces
machine-dependence and is off by default.

This is a **dense-matrix-only** method: it materialises and repeatedly
subsets the full \\n \times n\\ distance matrix, so it is suited to
instances small enough to hold that matrix. It offers no coordinate or
column-oracle path. For the matrix-free regime where the dense matrix is
infeasible, use
[`DropAddTSPoints()`](https://ms609.github.io/MaxMin/reference/DropAddTSPoints.md)
or [`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md)
(coordinate or distance-column oracle path), whose \\T_k\\ lands within
roughly a percent on the benchmark while scaling to far larger
instances.

## References

Resende MGC, Martí R, Gallego M, Duarte A (2010). “GRASP and path
relinking for the max-min diversity problem.” *Computers & Operations
Research*, **37**(3), 498–508.
[doi:10.1016/j.cor.2008.05.011](https://doi.org/10.1016/j.cor.2008.05.011)
.

## See also

[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
and
[`DropAddTSPoints()`](https://ms609.github.io/MaxMin/reference/DropAddTSPoints.md)
for scalable refinement;
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
for the proven optimum on small instances.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
res <- GraspPR(dist(pts), m = 5L, max_no_improve = 20L, elite_size = 4L,
               seed = 1L)
res$indices
#> [1]  3  4  5 24 25
```
