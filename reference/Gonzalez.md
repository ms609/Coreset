# Deterministic Gonzalez furthest-point selection

Greedy k-centre selection (González 1985) . Iteratively selects the
point furthest from the current selection, a 2-approximation to the
k-centre problem. The quality of the result depends on the first (seed)
point; by default `Gonzalez()` runs an **ensemble** of cheap `O(N)`
seeding strategies and keeps the selection with the largest minimum
pairwise distance
([`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)).
The default ensemble is three `"random_furthest"` starts – a
best-of-three selection. The random starts use the session RNG, so set a
seed ([`set.seed()`](https://rdrr.io/r/base/Random.html)) for a
reproducible selection. The deterministic `O(N)` anchors (`"centroid"`,
`"peripheral"`) and the costlier `O(N^2)` anchors (`"diameter"`,
`"anti_medoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`) are available as
opt-in `seed` strategies.

## Usage

``` r
Gonzalez(
  d = NULL,
  n,
  seed = .kDefaultEnsemble,
  pivots = NULL,
  points = NULL,
  N = NULL,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- d:

  A `dist` object, a square symmetric numeric matrix of pairwise
  distances, or a **distance-column oracle** function (see
  *Distance-column oracle*). Ignored when `points` is supplied.

- n:

  Integer: number of points to select. If `n >= N`, all indices are
  returned.

- seed:

  Integer or character (scalar or vector). An **integer** gives the
  explicit 1-based index of the first selected point (a single bare
  Gonzalez pass). A **length-1 character** names a single deterministic
  seeding strategy run as one bare pass: `"centroid"` (coordinates
  only), `"peripheral"` (two-sweep diameter-endpoint approximation),
  `"diameter"`, `"anti_medoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`, or
  `"first"` (index 1). A **length \> 1 character vector** – or the lone
  `"random_furthest"` token – requests an ensemble: each named anchor
  runs a full Gonzalez pass and the best result by
  [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) is
  returned with `strategy_results` and `winning_strategy` (character
  vector of all tied-best strategies) attributes. The
  `"random_furthest"` token expands to one start per element of
  `pivots`, labelled `random_furthest1`, `random_furthest2`, ...; named
  on its own it still runs the ensemble (one pass per pivot), so a
  single random start is best obtained via
  [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md).
  Valid ensemble anchors: any subset of
  `c("centroid", "peripheral", "random_furthest", "diameter", "anti_medoid", "medoid", "rowsum", "rownorm")`
  (`"centroid"` requires `points`). Default: `"random_furthest"` (three
  random starts; see `pivots`). See
  [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
  for anchor definitions. On the distance-column oracle path only an
  integer `seed` is honoured; a named or ensemble `seed` there warns and
  falls back to the peripheral seed (see *Distance-column oracle*).

- pivots:

  Integer vector of pivot indices over which the `"random_furthest"`
  ensemble token expands: each pivot contributes one start, seeded at
  the point furthest from it, so the vector's length sets the number of
  random-furthest starts. Left unspecified, three pivots are drawn with
  the session RNG (`sample.int(N, 3)`; set a seed for a reproducible
  selection). Pass `integer(0)`, `NA`, or `NULL` to disable the random
  starts, or an index vector to choose the pivots (and their count)
  explicitly. Disabling the random starts errors under the default
  `seed` (which names only `"random_furthest"`, leaving no anchor); pair
  it with a deterministic `seed` such as `"peripheral"`.

- points:

  Optional `N x dim` numeric coordinate matrix. When supplied, the
  selection is computed directly from coordinates in `O(N * n * dim)`
  time and `O(N)` memory, never materialising the `N x N` distance
  matrix (`d` is then unused). For Euclidean data the returned indices
  are identical to the matrix path. Only complete (non-`NA`) data is
  supported.

- N:

  Integer: the total number of elements. Required (and used) only on the
  distance-column oracle path, where it cannot be inferred from the
  closure; ignored for the matrix and coordinate paths.

- progress:

  Logical; show a progress bar during greedy selection on the
  distance-column oracle path (the only path slow enough to warrant
  one). Default: `TRUE` in interactive sessions, `FALSE` otherwise
  (`getOption("MaxMin.progress", interactive())`).

## Value

Integer vector of length `min(n, N)` of selected indices.

## Details

`Gonzalez()` accepts the distances in whichever of three forms suits the
data, all returning identical selections on the same metric:

- a **distance matrix** (`d`):

  a `dist` object or square matrix, held in full;

- a **coordinate matrix** (`points`):

  each needed distance is recomputed from coordinates on the fly in
  `O(N)` memory, never materialising the `N x N` matrix (Euclidean data
  only);

- a **distance-column oracle** (a function passed as `d`):

  for metrics with neither a stored matrix nor a coordinate embedding –
  e.g. tree-to-tree distances computed on demand. See *Distance-column
  oracle*.

## Distance-column oracle

When `d` is a function it is treated as a closure `colFn(i)` returning,
for a single 1-based index `i`, the length-`N` vector of distances from
element `i` to every element. Gonzalez calls it once per selected
element, maintaining a running nearest-distance vector, so the `N x N`
matrix is never built: `O(N * n)` oracle calls and `O(N)` memory. The
self-distance at position `i` may take any non-negative value; it is
masked before use. Because the count of elements cannot be inferred from
the closure, `N` must be supplied. Only an integer `seed` (a `first`
index) or a single deterministic two-sweep peripheral seed is reachable
from an oracle; the other anchors need either coordinates or `O(N^2)`
work. An explicitly named or ensemble `seed` on this path is ignored
with a warning and the peripheral seed is used instead.

## References

González TF (1985). “Clustering to minimize the maximum intercluster
distance.” *Theoretical Computer Science*, **38**, 293–306.
[doi:10.1016/0304-3975(85)90224-5](https://doi.org/10.1016/0304-3975%2885%2990224-5)
.

## See also

[`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
for the seed indices alone;
[`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md)
and
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
for higher-effort solvers.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
# Default: best of three random-furthest starts (set.seed for reproducibility):
Gonzalez(d, 5L)
#> [1] 14  4 26  5 28
#> attr(,"strategy_results")
#> attr(,"strategy_results")$random_furthest1
#> attr(,"strategy_results")$random_furthest1$s1
#> [1] 24
#> 
#> attr(,"strategy_results")$random_furthest1$idx
#> [1] 24  4  1  5 14
#> 
#> attr(,"strategy_results")$random_furthest1$t_k
#> [1] 1.701019
#> 
#> 
#> attr(,"strategy_results")$random_furthest2
#> attr(,"strategy_results")$random_furthest2$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$random_furthest2$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$random_furthest2$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest3
#> attr(,"strategy_results")$random_furthest3$s1
#> [1] 26
#> 
#> attr(,"strategy_results")$random_furthest3$idx
#> [1] 26 24 15  3 18
#> 
#> attr(,"strategy_results")$random_furthest3$t_k
#> [1] 1.468498
#> 
#> 
#> attr(,"winning_strategy")
#> [1] "random_furthest2"
# More random-furthest starts (length of `pivots` sets the count):
Gonzalez(d, 5L, pivots = sample.int(nrow(as.matrix(d)), 8L))
#> [1] 14  4 26  5 28
#> attr(,"strategy_results")
#> attr(,"strategy_results")$random_furthest1
#> attr(,"strategy_results")$random_furthest1$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$random_furthest1$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$random_furthest1$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest2
#> attr(,"strategy_results")$random_furthest2$s1
#> [1] 24
#> 
#> attr(,"strategy_results")$random_furthest2$idx
#> [1] 24  4  1  5 14
#> 
#> attr(,"strategy_results")$random_furthest2$t_k
#> [1] 1.701019
#> 
#> 
#> attr(,"strategy_results")$random_furthest3
#> attr(,"strategy_results")$random_furthest3$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$random_furthest3$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$random_furthest3$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest4
#> attr(,"strategy_results")$random_furthest4$s1
#> [1] 24
#> 
#> attr(,"strategy_results")$random_furthest4$idx
#> [1] 24  4  1  5 14
#> 
#> attr(,"strategy_results")$random_furthest4$t_k
#> [1] 1.701019
#> 
#> 
#> attr(,"strategy_results")$random_furthest5
#> attr(,"strategy_results")$random_furthest5$s1
#> [1] 5
#> 
#> attr(,"strategy_results")$random_furthest5$idx
#> [1]  5 26 14 21 28
#> 
#> attr(,"strategy_results")$random_furthest5$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest6
#> attr(,"strategy_results")$random_furthest6$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$random_furthest6$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$random_furthest6$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest7
#> attr(,"strategy_results")$random_furthest7$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$random_furthest7$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$random_furthest7$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest8
#> attr(,"strategy_results")$random_furthest8$s1
#> [1] 26
#> 
#> attr(,"strategy_results")$random_furthest8$idx
#> [1] 26 24 15  3 18
#> 
#> attr(,"strategy_results")$random_furthest8$t_k
#> [1] 1.468498
#> 
#> 
#> attr(,"winning_strategy")
#> [1] "random_furthest1" "random_furthest3" "random_furthest5" "random_furthest6"
#> [5] "random_furthest7"
# Or choose the pivots explicitly:
Gonzalez(d, 5L, pivots = c(1L, 10L, 20L))
#> [1]  5 26 14 21 28
#> attr(,"strategy_results")
#> attr(,"strategy_results")$random_furthest1
#> attr(,"strategy_results")$random_furthest1$s1
#> [1] 5
#> 
#> attr(,"strategy_results")$random_furthest1$idx
#> [1]  5 26 14 21 28
#> 
#> attr(,"strategy_results")$random_furthest1$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$random_furthest2
#> attr(,"strategy_results")$random_furthest2$s1
#> [1] 24
#> 
#> attr(,"strategy_results")$random_furthest2$idx
#> [1] 24  4  1  5 14
#> 
#> attr(,"strategy_results")$random_furthest2$t_k
#> [1] 1.701019
#> 
#> 
#> attr(,"strategy_results")$random_furthest3
#> attr(,"strategy_results")$random_furthest3$s1
#> [1] 24
#> 
#> attr(,"strategy_results")$random_furthest3$idx
#> [1] 24  4  1  5 14
#> 
#> attr(,"strategy_results")$random_furthest3$t_k
#> [1] 1.701019
#> 
#> 
#> attr(,"winning_strategy")
#> [1] "random_furthest1"
# Custom two-anchor ensemble:
Gonzalez(d, 5L, seed = c("diameter", "anti_medoid"))
#> [1] 14  4 26  5 28
#> attr(,"strategy_results")
#> attr(,"strategy_results")$diameter
#> attr(,"strategy_results")$diameter$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$diameter$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$diameter$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"strategy_results")$anti_medoid
#> attr(,"strategy_results")$anti_medoid$s1
#> [1] 14
#> 
#> attr(,"strategy_results")$anti_medoid$idx
#> [1] 14  4 26  5 28
#> 
#> attr(,"strategy_results")$anti_medoid$t_k
#> [1] 1.765223
#> 
#> 
#> attr(,"winning_strategy")
#> [1] "diameter"    "anti_medoid"
# A single strategy:
Gonzalez(d, 5L, seed = "diameter")
#> [1] 14  4 26  5 28
# An explicit start index (integer seed):
Gonzalez(d, 5L, seed = 1L)
#> [1]  1  5 24  4 14
# Matrix-free coordinate path (identical result, O(N) memory):
Gonzalez(n = 5L, points = pts, seed = 1L)
#> [1]  1  5 24  4 14
# Distance-column oracle: supply one column at a time, never the full matrix.
dm <- as.matrix(d)
colFn <- function(i) dm[, i]
identical(Gonzalez(colFn, 5L, N = nrow(dm), seed = 1L),
          Gonzalez(d, 5L, seed = 1L))
#> [1] TRUE
```
