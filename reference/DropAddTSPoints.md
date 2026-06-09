# Matrix-free DropAdd Tabu Search for the Max-Min Diversity Problem

Coordinate-based counterpart of
[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
that never materialises the dense \\n \times n\\ distance matrix. Each
needed distance column \\d(\cdot, x)\\ is recomputed from the supplied
coordinates on the fly in \\O(n \cdot \mathrm{dim})\\, giving \\O(n)\\
working memory. This lets the SOTA MaxMin heuristic of Porumbel et al.
(2011) run on point sets far larger than the matrix path can hold (R's
`as.matrix.dist` overflows at \\n = 46340\\; an \\n = 58000\\ double
matrix is roughly 27 GB).

## Usage

``` r
DropAddTSPoints(
  points,
  m,
  max_no_improve = 5000L,
  max_iter = NULL,
  time_budget_s = Inf,
  seed = NULL,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- points:

  A numeric \\n \times \mathrm{dim}\\ coordinate matrix (or an object
  coercible to one via `as.matrix`). Must be complete (no `NA`); the
  coordinate path is defined only for complete Euclidean data.

- m:

  Integer; subset size, \\2 \le m \le n\\.

- max_no_improve:

  Integer; stop after this many consecutive drop-add iterations that do
  not improve the best objective. The primary, deterministic stopping
  criterion. The search is RNG-free (ties broken by smallest index), so
  the result is reproducible and machine-independent. Default 5000.

- max_iter:

  Optional integer hard cap on main-loop iterations (excluding
  construction). `NULL` (default) leaves `max_no_improve` in sole
  control.

- time_budget_s:

  Optional wall-clock ceiling in seconds, checked periodically at
  iteration boundaries. Default `Inf` (no ceiling, fully reproducible).
  A finite value caps runtime but makes the result machine-dependent.

- seed:

  Optional integer; if non-`NULL`, `set.seed(seed)` is called at entry.
  The algorithm is deterministic up to ties (broken by smallest index),
  so the seed has no observable effect on the solution; it is exposed
  for API parity with stochastic methods and with
  [`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md).

- progress:

  Logical; show a start/done status line. Default: `TRUE` in interactive
  sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`).

## Value

List with elements

- indices:

  integer(m), 1-based selected indices sorted ascending.

- objective:

  numeric(1), achieved MaxMin objective \\\min\_{i \ne j \in S}
  d\_{ij}\\.

- secondary:

  numeric(1), achieved sum of pairwise distances over \\S\\
  (upper-triangle sum).

- time_s:

  numeric(1), wall-clock seconds spent.

- iters:

  integer(1), main-loop iterations executed (excluding the construction
  phase).

## Details

The construction (Algorithm 1), FIFO drop-add tabu search (Algorithm 2),
and streamlined neighbour evaluation (Algorithms 3-4) are identical to
[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md),
including the exclusion of the just-dropped point from the add
candidates for one iteration (Porumbel et al. (2011) , p.281). The MMDPo
objective optimised is \$\$\min\_{x,y \in X} d(x,y) + \epsilon
\sum\_{x,y \in X} d(x,y),\$\$ with \\\epsilon = 10^{-9}\\.

Seed deviation. Porumbel seeds the construction with the maximum-row-sum
point, \\\arg\max_x \sum_y d(x,y)\\, which costs \\O(n^2 \cdot
\mathrm{dim})\\ — the very expense this variant avoids. We substitute
the \\O(n \cdot \mathrm{dim})\\ proxy \\\arg\max_x \lVert x - \bar{x}
\rVert\\ (the point farthest from the coordinate centroid), which
approximates the peripheral max-row-sum point. The remainder of the
algorithm is faithful, so on instances where the two seed rules coincide
the trajectory matches
[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
exactly; otherwise the final quality is comparable (within a few percent
on the MaxMin objective, empirically).

Distances reproduce
[`stats::dist()`](https://rdrr.io/r/stats/dist.html)'s Euclidean bits
exactly (see `src/maximin_points.cpp` for the bit-equivalence argument),
so on data where the seeds coincide the matrix-free and matrix paths are
bit-identical under a toolchain whose `stats` package and this kernel
use the same floating-point contraction settings.

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.

## See also

[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
for the matrix-based path used on smaller instances.
