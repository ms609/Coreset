# DropAdd Tabu Search for the Max-Min Diversity Problem

Implements the DropAdd-TS algorithm of Porumbel et al. (2011) for
selecting a maximally-dispersed subset of `m` points from a distance
matrix. The procedure consists of a deterministic greedy construction
(Algorithm 1) followed by a FIFO drop-add tabu search (Algorithm 2) with
the streamlined neighbour-evaluation tricks of Algorithms 3 and 4.

## Usage

``` r
DropAddTS(
  d,
  m,
  max_no_improve = 5000L,
  max_iter = NULL,
  time_budget_s = Inf,
  seed = NULL,
  progress = getOption("MaxMin.progress", interactive()),
  .verify = FALSE,
  .trace = NULL
)
```

## Arguments

- d:

  A `dist` object or square symmetric numeric matrix.

- m:

  Integer; subset size, \\2 \le m \le n\\.

- max_no_improve:

  Integer; stop after this many consecutive drop-add iterations that do
  not improve the best objective. The primary, deterministic stopping
  criterion. The search is RNG-free (ties broken by smallest index), so
  for a given instance the result is reproducible and
  machine-independent. Default 5000.

- max_iter:

  Optional integer hard cap on iterations (excluding construction).
  `NULL` (default) leaves `max_no_improve` in sole control.

- time_budget_s:

  Optional wall-clock ceiling in seconds, checked at iteration
  boundaries. Default `Inf` (no ceiling, fully reproducible). A finite
  value caps runtime but makes the result machine-dependent.

- seed:

  Optional integer; if non-NULL, `set.seed(seed)` is called at entry.
  The algorithm is deterministic up to ties (broken by smallest index),
  so the seed has no observable effect on the solution; it is exposed
  for API parity with stochastic methods.

- progress:

  Logical; show a start/done status line. Default: `TRUE` in interactive
  sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`). No effect when
  `.verify = TRUE` (testing path).

- .verify:

  Logical (testing only); if `TRUE`, routes to the R reference loop and
  brute-force asserts the streamlined records at every iteration.
  Default `FALSE` (the C++ fast path).

- .trace:

  Optional environment (testing only); if supplied, the dropped and
  added index sequences are written into it as `drops` and `adds`.

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

The MMDPo objective optimised is \$\$\min\_{x,y \in X} d(x,y) + \epsilon
\sum\_{x,y \in X} d(x,y),\$\$ with \\\epsilon = 10^{-9}\\, so the
pairwise sum acts as a tie-break.

Tabu mechanics. The algorithm maintains an integer `iter` stamp for
every point, set to the iteration at which the point last entered or
left the selected set \\X\\. At each main-loop iteration the point with
the smallest `iter` value (the oldest member, FIFO) is dropped, and the
point in \\Z \setminus X\\ maximising lexicographically
\\(\mathrm{min\\dist}, \mathrm{sum\\dist})\\ is added. The FIFO
invariant guarantees that across any window of \\m\\ iterations every
initially-selected point is dropped exactly once before any re-eviction.

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.
