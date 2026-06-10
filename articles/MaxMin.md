# Introduction to MaxMin

The **Max-Min Diversity Problem** (MMDP) asks: given a fixed set of *N*
candidate items and a distance between every pair, choose *m* items so
that the closest pair in the selection is as far apart as possible. The
objective — often written T_(k) — is the minimum pairwise distance
within the chosen subset; a larger T_(k) means a more spread-out,
representative selection.

This problem arises naturally wherever you want a diverse sample from a
fixed pool: selecting field-survey sites to cover a landscape, picking
biological specimens for sequencing that span the available genetic
diversity, or choosing a representative subset of protein structures
from a database.

## Installation

``` r

install.packages("MaxMin")
```

## Quick start

`eurodist` is a built-in R `dist` object containing road distances (km)
between 21 European cities.

``` r

data(eurodist)

# Set a seed for a reproducible selection
set.seed(1)

# Select 4 maximally dispersed cities
idx <- FarFirst(eurodist, n = 4L)
MinDist(eurodist, idx)
#> [1] 2187

# View distances between chosen cities
as.matrix(eurodist)[idx, idx]
#>           Athens Lisbon Stockholm Milan
#> Athens         0   4532      3927  2282
#> Lisbon      4532      0      3231  2250
#> Stockholm   3927   3231         0  2187
#> Milan       2282   2250      2187     0
```

[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
returned the indices of the four cities whose nearest-neighbour distance
within the selection is largest.
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)
reports that value explicitly.

## Methods at a glance

MaxMin provides four solvers, and a minimum-distance calculator:

| Function | Quality | Speed | Stochastic? |
|----|----|----|----|
| [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) | Good (2-approximation) | Very fast | No |
| [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) | High (≈ 99 % optimal) | Fast | No |
| [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) | Highest | Moderate | Yes (`seed =`) |
| [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) | Optimal (NP-hard) | Slow | No |
| [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) | Scoring only | Instant | No |

We now walk through each of these in turn.

## Gonzalez: fast greedy selection

The González ([1985](#ref-Gonzalez1985)) algorithm builds the selection
greedily: start from a seed point, then repeatedly add whichever
unselected point is farthest from the current selection. This greedy
rule guarantees a 2-approximation to the optimal T_(k) and runs in O(*N*
· *m*) time.

By default
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
runs three `"random_furthest"` starts, each beginning from the point
furthest from a randomly selected ‘pivot’ point, and returns whichever
pass produced the highest T_(k).

``` r

set.seed(1)
picks <- FarFirst(eurodist, n = 6L)   # default: best of three random starts
```

More random starts, or pivots of your own devising, can be accomplished
by a vector of pivot points to a `pivots` vector.

``` r

set.seed(1)
nCities <- attr(eurodist, "Size")
# Use 5 random starts
pivots <- sample(nCities, 5)
picks <- FarFirst(eurodist, n = 6L, pivots = pivots)
```

Peripheral seeds may also be selected by deterministic methods: one or
more such methods can be selected via the `seed` argument. The best
solution found will be returned.

``` r

picks <- FarFirst(eurodist, n = 6L, seed = c("diameter", "anti_medoid"))
MinDist(eurodist, picks)
#> [1] 1014

attr(picks, "winning_strategy")
#> [1] "diameter"    "anti_medoid"
```

When only pairwise distances between objects are known — there are no
coordinates to work from — and computing or storing all N×N distances at
once would be too expensive, a function may be passed in place of the
distance matrix. This function will be passed one index `i`, and should
returns the distances from object `i` to all N objects — one column of
the matrix.
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
calls it $`n`$ times, so the full N×N matrix is never built. `N`, the
number of objects, must be supplied.

``` r

data("USArrests")
arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
StateDist <- function(i) {
  diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
  sqrt(rowSums(diffs ^ 2))
}
idx <- FarFirst(StateDist, n = 4L, N = nrow(arrestTypes), seed = 1L)
arrestTypes[idx, ]
#>                Murder Assault Rape
#> Alabama          13.2     236 21.2
#> North Dakota      0.8      45  7.3
#> North Carolina   13.0     337 16.1
#> Washington        4.0     145 26.2
```

## DropAdd tabu search

The **DropAdd** heuristic ([Porumbel et al., 2011](#ref-Porumbel2011))
refines an initial selection by alternately dropping and adding points
from the selection; it typically reaches ≈ 99 % of the optimal T_(k).

``` r

picksDA <- DropAdd(eurodist, m = 6L, plateau = 500L)

labels(eurodist)[picksDA]
#> [1] "Athens"    "Barcelona" "Cherbourg" "Lisbon"    "Stockholm" "Vienna"
attr(picksDA, "score")  # T_k achieved
#> [1] 1294
attr(picksDA, "iters")      # iterations completed
#> [1] 516
attr(picksDA, "time_s")     # wall-clock seconds
#> [1] 0.001
```

The algorithm terminates after `plateau` iterations do not improve
T_(k). Where *N* is too large for a distance matrix to fit in memory
(roughly N \> 46 000), pass a coordinate matrix via
`DropAdd(points = ...)`.

## GRASP with path relinking

**GRASP + path relinking** ([Resende et al., 2010](#ref-Resende2010))
combines a *randomised* greedy construction phase with extended local
search and then refines an elite set of good solutions by interpolating
between elite-pair trajectories (*path relinking*). It achieves the
highest T_(k) of the three heuristics, at a proportionally higher cost.
Because the construction phase is stochastic, always set `seed` for
reproducibility:

``` r

res_gr <- Grasp(eurodist, m = 6L, plateau = 50L, seed = 42L)

labels(eurodist)[res_gr]
#> [1] "Athens"    "Barcelona" "Calais"    "Lisbon"    "Stockholm" "Vienna"
attr(res_gr, "score")
#> [1] 1305
attr(res_gr, "pr_calls")   # path-relinking calls performed
#> [1] 90
```

`plateau` controls how many consecutive non-improving GRASP iterations
trigger termination; `timeBudgetS` is available as an absolute time cap
for interactive or batch use.

## Comparing methods on a simulated example

To see how the methods relate visually, we generate 50 points in two
dimensions and select *m* = 8 from each.

``` r

set.seed(42)
pts <- matrix(rnorm(100), ncol = 2)   # 50 points, 2 dimensions
d50 <- dist(pts)
m   <- 8L
```

``` r

idx_ff <- FarFirst(d50, n = m)
res_da50 <- DropAdd(d50, m = m, plateau = 500L)
res_gr50 <- Grasp(d50, m = m, plateau = 50L, seed = 42L)
```

Even a small difference in T_(k) can correspond to a meaningfully more
dispersed selection.

``` r

scores <- c(
  FarFirst  = MinDist(d50, idx_ff),
  DropAdd = attr(res_da50, "score"),
  Grasp   = attr(res_gr50, "score")
)
round(scores, 3)
#> FarFirst  DropAdd    Grasp 
#>    1.416    1.428    1.428
```

Plotting the selections against the point cloud makes the differences
concrete. When two methods select the same index, their symbols overlap;
the T_(k) table above captures the quality distinction even when the
visual overlap is high.

``` r

methods <- list(
  FarFirst  = idx_ff,
  DropAdd = res_da50,
  Grasp   = res_gr50
)
cols <- c(FarFirst = "#E41A1C", DropAdd = "#377EB8", Grasp = "#4DAF4A")
pchs <- c(FarFirst = 24L, DropAdd = 21L, Grasp = 22L)
# Small per-method shift so coincident selections remain distinguishable
dx   <- c(FarFirst = 0, DropAdd =  0.04, Grasp = -0.04)
dy   <- c(FarFirst = 0, DropAdd = -0.04, Grasp =  0.04)

plot(pts, pch = 1L, col = "grey75", asp = 1L,
     xlab = "x", ylab = "y",
     main = "MaxMin method comparison")

for (nm in names(methods)) {
  sel <- methods[[nm]]
  points(
    pts[sel, 1L] + dx[nm],
    pts[sel, 2L] + dy[nm],
    pch = pchs[nm], col = cols[nm], bg = cols[nm], cex = 1.6
  )
}

legend("topright",
       legend = paste0(names(methods), "  (Tₖ = ", round(scores, 3), ")"),
       pch = pchs, col = cols, pt.bg = cols, pt.cex = 1.4, bty = "n")
```

![](MaxMin_files/figure-html/compare-plot-1.png)

Selections returned by each method on 50 random 2-D points (m = 8).
Coloured symbols mark selected points; grey circles are the full
candidate set. A small jitter separates symbols at shared indices.

## Exact solution

For small instances (roughly N ≤ 25–30),
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
solves the problem to proven optimality via a node-packing integer
programme ([Sayyady & Fathi, 2016](#ref-Sayyady2016)), using the
**highs** solver.

First we must install **highs**:

``` r

install.packages("highs")
```

``` r

set.seed(1L)
pts30 <- matrix(rnorm(60L), ncol = 2L)
d30   <- dist(pts30)
```

``` r

res_ex <- ExactMaxMin(d30, m = 6L, timeBudgetS = 30L)

res_ex$proven      # TRUE  ⟹  objective is the global optimum
#> [1] TRUE
res_ex$objective
#> [1] 1.616311

# Compare to the greedy heuristic on the same instance
res_ff30 <- FarFirst(d30, n = 6L)
c(exact    = res_ex$objective,
  Gonzalez = MinDist(d30, res_ff30))
#>    exact Gonzalez 
#> 1.616311 1.616311
```

`$proven = TRUE` certifies that no selection can achieve a higher T_(k).
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
is NP-hard; instances with more than ~30 candidates require a time
budget or should be tackled with the heuristics above.

## Scoring

### `MinDist()`

[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)
computes the T_(k) objective for any index set. It accepts a `dist`
object, a square distance matrix, or a coordinate matrix via the
`points` argument:

``` r

MinDist(d50, idx_ff)                               # from dist
#> [1] 1.416314
MinDist(as.matrix(d50), idx_ff)                    # from square matrix
#> [1] 1.416314
MinDist(points = pts, idx = idx_ff)                # from coordinates
#> [1] 1.416314
```

## When to use which method

| Scenario | Recommended |
|----|----|
| Speed matters most, N up to a few thousand | [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) (ensemble default) |
| Deterministic, reproducible refinement | [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) |
| Best quality, set `seed` for reproducibility | `Grasp(seed = ...)` |
| N \> 46 000 (distance matrix infeasible) | `DropAdd(points = ...)` or `FarFirst(points = ...)` |
| Arbitrary metric with no coordinate embedding | `FarFirst(<column function>, N = ...)` |
| Proven optimum, N ≤ ~ 25–30, **highs** installed | [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) |

The heuristics are complementary:
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) is
O(*N* · *m*) and deterministic — an instant first result.
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) is
deterministic and reproducible (no RNG), typically reaching ≈ 99 % of
optimal; the `seed` parameter governs only the RNG used for
tie-breaking, not the algorithm itself.
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) is
stochastic and usually edges out
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) on
T_(k), but requires a `seed` for reproducibility.

## Related packages

The CRAN package
[**maximin**](https://CRAN.R-project.org/package=maximin) constructs
continuous space-filling designs by generating new points.

## References

González, T. F. (1985). Clustering to minimize the maximum intercluster
distance. *Theoretical Computer Science*, *38*, 293–306.
<https://doi.org/10.1016/0304-3975(85)90224-5>

Porumbel, D., Hao, J.-K., & Glover, F. (2011). A simple and effective
algorithm for the MaxMin diversity problem. *Annals of Operations
Research*, *186*, 275–293. <https://doi.org/10.1007/s10479-011-0898-z>

Resende, M. G. C., Martí, R., Gallego, M., & Duarte, A. (2010). GRASP
and path relinking for the max-min diversity problem. *Computers &
Operations Research*, *37*(3), 498–508.
<https://doi.org/10.1016/j.cor.2008.05.011>

Sayyady, F., & Fathi, Y. (2016). An integer programming approach for
solving the p-dispersion problem. *European Journal of Operational
Research*, *253*(1), 216–225.
<https://doi.org/10.1016/j.ejor.2016.02.026>
