# Greedy farthest-first point selection

Greedy farthest-first selection (González 1985; Hochbaum and Shmoys
1985) iteratively selects the point furthest from the current selection
to yield a 2-approximation to the *k*-centre and Max Min Diversity
problems.

## Usage

``` r
FarFirst(
  k,
  d = NULL,
  points = NULL,
  N = NULL,
  strategy = "random_furthest",
  nSeeds = 3L
)
```

## Arguments

- k:

  Integer: number of points to select.

- d:

  A `dist` object, a square numeric matrix of pairwise distances, or a
  distance function that takes an index `i` and returns the distance
  from `i` to each other element (optionally including the
  self-distance). Ignored when `points` is supplied.

- points:

  Optional `N x dim` numeric coordinate matrix. When supplied, the
  selection is computed directly from coordinates in \\O(N \cdot k \cdot
  dim)\\ time and \\O(N)\\ memory, which avoids creating an \\O(N^2)\\
  distance matrix. Missing entries (`NA`) are not supported.

- N:

  Integer: the total number of elements. Required (and used) only on the
  distance-column oracle path, where it cannot be inferred from the
  closure; ignored for the matrix and coordinate paths.

- strategy:

  Integer or character defining how to seed the greedy pass. Pass the
  name of one or more seeding strategies described in
  [`PickPoint()`](https://ms609.github.io/MaxMin/reference/PickPoint.md)
  to run each strategy and return the best solution.

- nSeeds:

  Integer: number of distinct seeds to draw under the (default)
  `"random_furthest"` strategy. Three starts captures most of the gain
  from restarting — the improvement curve bends early (knee at n ≈ 3–4
  across benchmarks) and additional restarts add little. For
  higher-quality solutions, prefer
  [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md):
  tabu search escapes the farthest-first construction family where
  restarts plateau.

## Value

`FarFirst()` returns an integer vector with class `MaxMinSelection`,
listing the selected indices in the order they were selected. Attributes
report:

- `score`: the selection's minimum pairwise distance (\\T_k\\).

- `winning_strategy`: character vector listing strategies that attained
  the optimal score.

- `strategy_results`: results for each strategy.

## Progress bar

In interactive sessions, the distance-column path shows a progress bar.
To toggle, set `options("MaxMin.progress" = FALSE)` (or `TRUE`).

## Parallelism

To parallelize computation when OpenMP is available, set the
`"mc.cores"` option:

    options(mc.cores = 2L)                       # use a fixed number of cores
    options(mc.cores = parallel::detectCores())  # or all available cores

The main use case for parallelization is when `nSeeds` is a multiple of
`mc.cores`, so each seed point can be evaluated in parallel.

## References

González TF (1985). “Clustering to minimize the maximum intercluster
distance.” *Theoretical Computer Science*, **38**, 293–306.
[doi:10.1016/0304-3975(85)90224-5](https://doi.org/10.1016/0304-3975%2885%2990224-5)
.  
  
Hochbaum DS, Shmoys DB (1985). “A best possible heuristic for the
\\k\\-center problem.” *Mathematics of Operations Research*, **10**(2),
180–184.
[doi:10.1287/moor.10.2.180](https://doi.org/10.1287/moor.10.2.180) .

## See also

[`PickPoint()`](https://ms609.github.io/MaxMin/reference/PickPoint.md)
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
FarFirst(5L, d)
#> 5 elements (14 4 26 5 28) selected by farthest-first (best of 3 strategies, winner random_furthest1), each at distance >= 1.765

# More random-furthest starts:
FarFirst(5L, d, nSeeds = 15L)
#> 5 elements (4 14 26 5 28) selected by farthest-first (best of 5 strategies, 3 tied: random_furthest1, random_furthest2, random_furthest3), each at distance >= 1.765

# Custom two-anchor ensemble:
FarFirst(5L, d, strategy = c("diameter", "anti_medoid"))
#> 5 elements (14 4 26 5 28) selected by farthest-first (best of 2 strategies, 2 tied: diameter, anti_medoid), each at distance >= 1.765

# A single strategy:
FarFirst(5L, d, strategy = "diameter")
#> 5 elements (14 4 26 5 28) selected by farthest-first, each at distance >= 1.765

# An explicit start index (integer strategy):
FarFirst(5L, d, strategy = 1L)
#> 5 elements (1 5 24 4 14) selected by farthest-first, each at distance >= 1.701

# Matrix-free coordinate path (identical result, O(N) memory):
FarFirst(5L, points = pts, strategy = 1L)
#> 5 elements (1 5 24 4 14) selected by farthest-first, each at distance >= 1.701

# Distance-column oracle: supply one column at a time, never the full matrix.
data("USArrests")
arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
StateDist <- function(i) {
  diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
  sqrt(rowSums(diffs ^ 2))
}
idx <- FarFirst(4L, StateDist, N = nrow(arrestTypes), strategy = 1L)
arrestTypes[idx, ]
#>                Murder Assault Rape
#> Alabama          13.2     236 21.2
#> North Dakota      0.8      45  7.3
#> North Carolina   13.0     337 16.1
#> Washington        4.0     145 26.2
```
