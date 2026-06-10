# Deterministic Gonzalez furthest-point selection

Greedy k-centre selection (González 1985) . Iteratively selects the
point furthest from the current selection, a 2-approximation to the
k-centre problem. The quality of the result depends on the first (seed)
point; by default `FarFirst()` runs three starts from randomly selected
peripheral seeds. The deterministic `O(N)` anchors (`"centroid"`,
`"peripheral"`) and the costlier `O(N^2)` anchors (`"diameter"`,
`"anti_medoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`) are alternative
`seed` strategies.

## Usage

``` r
FarFirst(
  d = NULL,
  m,
  method = .kDefaultEnsemble,
  pivots = NULL,
  points = NULL,
  N = NULL,
  progress = getOption("MaxMin.progress", interactive())
)
```

## Arguments

- d:

  A `dist` object, a square numeric matrix of pairwise distances, or a
  distance function (see *§Distance function*). Asymmetric matrices are
  accepted; symmetry is not checked (an `O(N^2)` check is intentionally
  omitted), and the algorithm treats \\d\_{ij}\\ and \\d\_{ji}\\ as
  independent. Ignored when `points` is supplied.

- m:

  Integer: number of points to select. If `m > N`, all `N` indices are
  returned in Gonzalez (farthest-first) order.

- method:

  Integer or character (scalar or vector); how to seed the greedy pass
  (matching the `method` argument of
  [`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)).
  An **integer** gives the explicit 1-based index of the first selected
  point (a single bare Gonzalez pass). A **length-1 character** names a
  single deterministic seeding strategy run as one bare pass:
  `"centroid"` (coordinates only), `"peripheral"` (two-sweep
  diameter-endpoint approximation), `"diameter"`, `"anti_medoid"`,
  `"medoid"`, `"rowsum"`, `"rownorm"`, or `"first"` (index 1). A
  **length \> 1 character vector** – or the lone `"random_furthest"`
  token – requests an ensemble: each named anchor runs a full Gonzalez
  pass and the best result by
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
  integer `method` is honoured; a named or ensemble `method` there warns
  and falls back to the peripheral seed (see *Distance-column oracle*).

- pivots:

  Integer vector of pivot indices over which the `"random_furthest"`
  ensemble token expands: each pivot contributes one start, seeded at
  the point furthest from it, so the vector's length sets the number of
  random-furthest starts. Left unspecified, three pivots are drawn with
  the session RNG (`sample.int(N, 3)`; set a seed for a reproducible
  selection). Pass `integer(0)`, `NA`, or `NULL` to disable the random
  starts, or an index vector to choose the pivots (and their count)
  explicitly. Disabling the random starts errors under the default
  `method` (which names only `"random_furthest"`, leaving no anchor);
  pair it with a deterministic `method` such as `"peripheral"`.

- points:

  Optional `N x dim` numeric coordinate matrix. When supplied, the
  selection is computed directly from coordinates in `O(N * m * dim)`
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

Integer vector of length `min(m, N)` of selected indices, in
farthest-first (greedy) selection order – **not** sorted (unlike
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md),
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) and
`ExactMaxMin()$indices`, which return ascending indices). The achieved
\\T_k\\ (the selection's minimum pairwise distance) is attached as
attribute `score`. An ensemble `method` additionally carries
`strategy_results` and `winning_strategy` attributes.

## Details

Distances may be provided as:

- a **distance matrix** (`d`):

  a `dist` object or square matrix, held in full;

- a **coordinate matrix** (`points`):

  each needed distance is recomputed from coordinates on the fly in
  `O(N)` memory (Euclidean data only);

- a distance function (function passed as `d`):

  for metrics with neither a stored matrix nor a coordinate embedding,
  where distances may be computed on demand.

## Distance function

When `d` is a function, it will be passed a single 1-based index `i`,
and should return the distances from element `i` to every element in
turn, optionally omitting entry `i`, the self-distance. The function
will be called once per selected element, to avoid building a complete
`N x N` matrix.

## References

González TF (1985). “Clustering to minimize the maximum intercluster
distance.” *Theoretical Computer Science*, **38**, 293–306.
[doi:10.1016/0304-3975(85)90224-5](https://doi.org/10.1016/0304-3975%2885%2990224-5)
.

## See also

[`MaxMinSeed()`](https://ms609.github.io/MaxMin/reference/MaxMinSeed.md)
for the seed indices alone;
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) and
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
for higher-effort solvers.

## Examples

``` r
set.seed(1)
pts <- matrix(rnorm(60), ncol = 2)
d <- dist(pts)
# Default: best of three random-furthest starts (set.seed for reproducibility):
FarFirst(d, 5L)
#> [1] 14  4 26  5 28
#> attr(,"score")
#> [1] 1.765223
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
FarFirst(d, 5L, pivots = sample.int(nrow(as.matrix(d)), 8L))
#> [1] 14  4 26  5 28
#> attr(,"score")
#> [1] 1.765223
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
FarFirst(d, 5L, pivots = c(1L, 10L, 20L))
#> [1]  5 26 14 21 28
#> attr(,"score")
#> [1] 1.765223
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
FarFirst(d, 5L, method = c("diameter", "anti_medoid"))
#> [1] 14  4 26  5 28
#> attr(,"score")
#> [1] 1.765223
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
FarFirst(d, 5L, method = "diameter")
#> [1] 14  4 26  5 28
#> attr(,"score")
#> [1] 1.765223
# An explicit start index (integer method):
FarFirst(d, 5L, method = 1L)
#> [1]  1  5 24  4 14
#> attr(,"score")
#> [1] 1.701019
# Matrix-free coordinate path (identical result, O(N) memory):
FarFirst(m = 5L, points = pts, method = 1L)
#> [1]  1  5 24  4 14
#> attr(,"score")
#> [1] 1.701019

# Distance-column oracle: supply one column at a time, never the full matrix.
data("USArrests")
arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
StateDist <- function(i) {
  diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
  sqrt(rowSums(diffs ^ 2))
}
idx <- FarFirst(StateDist, m = 4L, N = nrow(arrestTypes), method = 1L)
arrestTypes[idx, ]
#>                Murder Assault Rape
#> Alabama          13.2     236 21.2
#> North Dakota      0.8      45  7.3
#> North Carolina   13.0     337 16.1
#> Washington        4.0     145 26.2
```
