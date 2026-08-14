# DropAdd Tabu Search for the Max-Min Diversity Problem

`DropAdd()` selects a maximally-dispersed subset of `k` points using the
DropAdd tabu search algorithm, which comprises a greedy construction
followed by a first-in, first-out drop-add tabu search, with streamlined
neighbour-evaluation tricks (algorithms 1–4 in Porumbel et al. 2011) .

## Usage

``` r
DropAdd(
  k,
  d = NULL,
  plateau = 5000L,
  maxSeconds = Inf,
  points = NULL,
  maxCandidates = 46340L,
  seed = NULL,
  N = NULL
)
```

## Arguments

- k:

  Integer: subset size, \\2 \le k \le N\\.

- d:

  A `dist` object, a square symmetric numeric matrix, or a
  distance-column function (see §*Distance-column function*).

- plateau:

  Integer: stop after this many consecutive drop-add iterations do not
  improve the score.

- maxSeconds:

  Numeric: terminate search after this many seconds have elapsed.

- points:

  A numeric \\N \times \mathrm{dim}\\ coordinate matrix (or an object
  coercible to one via `as.matrix`). Must be complete (no `NA`). Ignored
  if `d` specified. Avoids creating an \\N \times N\\ distance matrix,
  enabling use at \\N \ge 46340\\).

- maxCandidates:

  Integer: a composable-coreset tractability cap. When the number of
  candidate points `N` exceeds `maxCandidates`,
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
  thins the candidates to a `maxCandidates`-point coreset (with the
  deterministic, RNG-free `"peripheral"` seed, so no random stream is
  perturbed), the solver runs on the coreset, and the chosen indices are
  mapped back to the original numbering. This lets the solver produce a
  solution at scales where it would otherwise be intractable.
  `maxCandidates = 0` (or `Inf`) disables thinning and runs on the full
  problem; a cap at or above `N` is a no-op. A cap below `k` is an
  error. The default is `46340L`, the dense-distance-matrix feasibility
  ceiling Thinning is **on by default**: an input larger than the cap is
  thinned (and a warning is emitted) unless `maxCandidates = 0` is
  passed.

- seed:

  Optional integer: a 1-based start index that overrides the
  construction's default warm-start seed. `NULL` (default) keeps the
  method's own seed. Not supported when `maxCandidates = 0L`.

- N:

  Integer: the total number of elements. Required only if `d` is a
  function.

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
[`print.MaxMinSelection()`](https://ms609.github.io/MaxMin/reference/print.MaxMin.md));
it is otherwise an ordinary integer vector.

## Progress bar

In interactive sessions, status messages are shown. To toggle, set
`options("MaxMin.progress" = FALSE)` (or `TRUE`).

## Parallelism

To parallelize computation when OpenMP is available, set the
`"mc.cores"` option:

    options(mc.cores = 2L)                       # use a fixed number of cores
    options(mc.cores = parallel::detectCores())  # or all available cores

## Distance function

When `d` is a function, `d(i)` must return the distances from element
`i` to every element (length `N`, with the self-distance ignored) or to
every *other* element (length `N - 1`, in order). `N` is required, and
memory is \\O(N)\\. This suits metrics where no stored matrix or
coordinate embedding is available.

It is likely that `d` will be called many times; unless `d` implements
caching, specifying a distance matrix is likely to require less
calculation than the multiple calls to `d`, where memory permits.

## References

Porumbel D, Hao J, Glover F (2011). “A simple and effective algorithm
for the MaxMin diversity problem.” *Annals of Operations Research*,
**186**, 275–293.
[doi:10.1007/s10479-011-0898-z](https://doi.org/10.1007/s10479-011-0898-z)
.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(200), ncol = 2)
DropAdd(5L, dist(pts))
#> 5 elements (14 41 56 87 97) selected by DropAdd tabu search, each at distance >= 2.293

# Composable coreset: thin to 40 candidates with farthest-first, then run
# DropAdd on the coreset. Returned indices are original-space row indices.
suppressWarnings(DropAdd(5L, points = pts, maxCandidates = 40L))
#> 5 elements (14 41 56 87 97) selected by DropAdd tabu search, each at distance >= 2.293

# Disable thinning on the full problem
DropAdd(5L, points = pts, maxCandidates = 0L)
#> 5 elements (14 41 56 87 97) selected by DropAdd tabu search, each at distance >= 2.293

# Distance function; `cache` memoizes the columns to reduce computation.
data("USArrests")
ArrestDist <- function(dat) {
  scaled <- scale(as.matrix(dat))                 # derived once
  cache <- new.env(parent = emptyenv())
  function(i) {
    key <- as.character(i)
    if (is.null(cache[[key]])) {
      cache[[key]] <- sqrt(rowSums(sweep(scaled, 2, scaled[i, ], "-") ^ 2))
    }
    cache[[key]]
  }
}
arrests <- USArrests[, c("Murder", "Assault", "Rape")]
idx <- DropAdd(4L, ArrestDist(arrests), N = nrow(arrests), plateau = 200L)
USArrests[idx, ]
#>              Murder Assault UrbanPop Rape
#> Alaska         10.0     263       48 44.5
#> Delaware        5.9     238       72 15.8
#> Georgia        17.4     211       60 25.8
#> North Dakota    0.8      45       44  7.3
```
