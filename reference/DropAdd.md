# DropAdd Tabu Search for the Max-Min Diversity Problem

Implements the DropAdd-TS algorithm of Porumbel et al. (2011) for
selecting a maximally-dispersed subset of `m` points from a distance
matrix. The procedure consists of a deterministic greedy construction
(Algorithm 1) followed by a FIFO drop-add tabu search (Algorithm 2) with
the streamlined neighbour-evaluation tricks of Algorithms 3 and 4.

## Usage

``` r
DropAdd(
  d = NULL,
  k,
  plateau = 5000L,
  maxIter = NULL,
  maxSeconds = Inf,
  progress = getOption("MaxMin.progress", interactive()),
  points = NULL,
  .verify = FALSE,
  .trace = NULL
)
```

## Arguments

- d:

  A `dist` object or square symmetric numeric matrix. Mutually exclusive
  with `points`; supply exactly one.

- k:

  Integer; subset size, \\2 \le k \le n\\.

- plateau:

  Integer; stop after this many consecutive drop-add iterations that do
  not improve the best objective. The primary, deterministic stopping
  criterion. The search is RNG-free (ties broken by smallest index), so
  for a given instance the result is reproducible and
  machine-independent. Default 5000.

- maxIter:

  Optional integer hard cap on iterations (excluding construction).
  `NULL` (default) leaves `plateau` in sole control.

- maxSeconds:

  Optional wall-clock ceiling in seconds, checked at iteration
  boundaries. Default `Inf` (no ceiling, fully reproducible). A finite
  value caps runtime but makes the result machine-dependent.

- progress:

  Logical; show a start/done status line. Default: `TRUE` in interactive
  sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`). No effect when
  `.verify = TRUE` (testing path).

- points:

  A numeric \\n \times \mathrm{dim}\\ coordinate matrix (or an object
  coercible to one via `as.matrix`). Must be complete (no `NA`).
  Mutually exclusive with `d`; supply exactly one. When supplied, the
  algorithm never materialises the dense \\n \times n\\ distance matrix,
  giving \\O(n)\\ working memory and enabling use at \\n\\ far exceeding
  the matrix path's ceiling (R's `as.matrix.dist` overflows at \\n =
  46340\\).

- .verify:

  Logical (testing only); if `TRUE`, routes to the R reference loop and
  brute-force asserts the streamlined records at every iteration.
  Default `FALSE` (the C++ fast path). Only applies to the `d` path.

- .trace:

  Optional environment (testing only); if supplied, the dropped and
  added index sequences are written into it as `drops` and `adds`.

## Value

`DropAdd()` returns an integer vector of length `k` containing the
1-based selected indices **sorted ascending** (unlike
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md),
which returns farthest-first order), with attributes:

- score:

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

The vector has class `"MaxMinSelection"` and prints as a one-line
summary (see
[print.MaxMinSelection](https://ms609.github.io/MaxMin/reference/print.MaxMin.md));
it is otherwise an ordinary integer vector.

## Details

The MMDPo objective optimised is \$\$\min\_{x,y \in X} d(x,y) + \epsilon
\sum\_{x,y \in X} d(x,y),\$\$ with \\\epsilon = 10^{-9}\\, so the
pairwise sum acts as a tie-break.

Tabu mechanics. The algorithm maintains an integer `iter` stamp for
every point, set to the iteration at which the point last entered or
left the selected set \\X\\. At each main-loop iteration the point with
the smallest `iter` value (the oldest member, FIFO) is dropped, and the
point in \\Z \setminus X\\ maximising lexicographically
\\(\mathrm{minDist}, \mathrm{sumDist})\\ is added. The FIFO invariant
guarantees that across any window of \\m\\ iterations every
initially-selected point is dropped exactly once before any re-eviction.

Time budget behaviour. The time budget (`maxSeconds`) is checked at most
once every 256 iterations (matrix-free path) or 1024 iterations (matrix
path). On large instances where each iteration is slow, the actual
elapsed time may exceed the specified budget by up to one iteration's
worth of computation.

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.
