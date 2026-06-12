# Introduction to MaxMin

The MaxMin package selects a subset that represents a fixed pool of *N*
items, based on one of two complementary objectives:

The **Max-Min Diversity Problem** (MMDP, the discrete *p*-dispersion
objective) selects $`k`$ elements such that the minimum distance between
any pair of selected elements is as large as possible; the chosen
elements are maximally separated. This can reward selections that leave
the interior of the set unrepresented.

This objective is suited to defining a representative sample from a
fixed pool: picking biological specimens for sequencing that span
available diversity, or choosing a representative subset of protein
structures from a database.

The **discrete *k*-centre problem** selects $`k`$ elements such that the
maximum distance from any element in the original set to a selected
element is as small as possible. In ensuring that each point has a
nearby representative, this objective can select points that reflect a
central compromise, rather than selections that are closer to more local
points.

This objective is useful when selecting centres that represent each
point in a dataset: for example, siting fire stations to guarantee that
all buildings can be reached within a given response time.

An approximate selection that satisfies both objectives within a factor
of two of their respective optima can be attained by a greedy
farthest-first algorithm ([González, 1985](#ref-Gonzalez1985)), though
the exact optima typically differ; dispersion spreads to the extremes,
whereas covering reaches into the interior. MaxMin provides approximate
and exact solvers for each objective.

## Installation

``` r

install.packages("MaxMin")
```

## Quick start

Our examples employ a built-in R `dist` object that contains road
distances (in km) between 21 European cities.

``` r

# Load the `eurodist` dist object
data("eurodist")

# Set a seed for a reproducible selection
set.seed(1)

# Select 4 maximally dispersed cities
ffPick <- FarFirst(eurodist, k = 4L)

# View distances between selected cities
as.matrix(eurodist)[ffPick, ffPick]
#>           Athens Lisbon Stockholm Milan
#> Athens         0   4532      3927  2282
#> Lisbon      4532      0      3231  2250
#> Stockholm   3927   3231         0  2187
#> Milan       2282   2250      2187     0

# Quickly extract the minimum distance for a given selection
MinDist(eurodist, ffPick)
#> [1] 2187
```

[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
returned the indices of the four cities whose nearest-neighbour distance
within the selection is largest.
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)
reports that value explicitly.

## Methods at a glance

MaxMin provides solvers for each objective:

| Function | Objective | Quality | Speed | Stochastic? |
|----|----|----|----|----|
| [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) | both | Good (2-approximation) | Very fast | No |
| [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) | MMDP | High (≈ 99 % optimal) | Fast | No |
| [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) | MMDP | Highest | Moderate | Yes ([`set.seed()`](https://rdrr.io/r/base/Random.html)) |
| [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) | MMDP | Optimal (NP-hard) | Slow | No |
| [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) | k-centre | Near-optimal (CDSh) | Fast | No |
| [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md) | k-centre | Optimal (NP-hard) | Slow | No |

To compute the score for an arbitrary selection of points under each
objective, use
[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)
(MMDP) and
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
(k-centre).

## Fast greedy selection

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
picks <- FarFirst(eurodist, k = 6L)   # default: best of eight random starts
```

More random starts can be requested via the `nseeds` argument.

``` r

set.seed(1)
picks <- FarFirst(eurodist, k = 6L, nseeds = 12L)
```

Peripheral seeds may also be selected by deterministic methods: one or
more such methods can be selected via the `method` argument. The best
solution found will be returned.

``` r

picks <- FarFirst(eurodist, k = 6L, method = c("diameter", "anti_medoid"))
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
calls it $`k`$ times, so the full N×N matrix is never built. `N`, the
number of objects, must be supplied.

``` r

data("USArrests")
arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
StateDist <- function(i) {
  diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
  sqrt(rowSums(diffs ^ 2))
}
idx <- FarFirst(StateDist, k = 4L, N = nrow(arrestTypes), method = 1L)
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

picksDA <- DropAdd(eurodist, k = 6L, plateau = 500L)

picksDA
#> 6 elements (1 2 5 12 20 21) selected by DropAdd tabu search, each at distance >= 1294
labels(eurodist)[picksDA]
#> [1] "Athens"    "Barcelona" "Cherbourg" "Lisbon"    "Stockholm" "Vienna"
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
Because the construction phase is stochastic, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) for a
reproducible run:

``` r

set.seed(42)
res_gr <- Grasp(eurodist, k = 6L, plateau = 50L)
res_gr
#> 6 elements (1 2 4 12 20 21) selected by GRASP with path-relinking, each at distance >= 1305
labels(eurodist)[res_gr]
#> [1] "Athens"    "Barcelona" "Calais"    "Lisbon"    "Stockholm" "Vienna"
attr(res_gr, "pr_calls")   # path-relinking calls performed
#> [1] 90
```

`plateau` controls how many consecutive non-improving GRASP iterations
trigger termination; `maxSeconds` is available as an absolute time cap
for interactive or batch use.

## Comparing methods on a simulated example

To see how the methods relate visually, we generate 50 points in two
dimensions and select *k* = 8 from each.

``` r

set.seed(42)
pts <- matrix(rnorm(100), ncol = 2)   # 50 points, 2 dimensions
d50 <- dist(pts)
k   <- 8L
```

``` r

idx_ff <- FarFirst(d50, k)
res_da50 <- DropAdd(d50, k = k, plateau = 500L)
set.seed(42)
res_gr50 <- Grasp(d50, k = k, plateau = 50L)
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

Selections returned by each method on 50 random 2-D points (k = 8).
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

res_ex <- ExactMaxMin(d30, k = 6L, maxSeconds = 30L)

res_ex$proven      # TRUE  ⟹  objective is the global optimum
#> [1] TRUE
res_ex$objective
#> [1] 1.616311

# Compare to the greedy heuristic on the same instance
res_ff30 <- FarFirst(d30, k = 6L)
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

## Covering: the k-centre problem

Every method above pursues *dispersion* — it spreads the selection so
its members are mutually far apart. The **k-centre** problem instead
pursues *coverage*: it minimises the covering radius *R*, the largest
distance from any point to its nearest chosen centre, so that no point
of the pool is left far from a representative ([González,
1985](#ref-Gonzalez1985); [Hochbaum & Shmoys, 1985](#ref-Hochbaum1985)).
Where dispersion reaches for the extremes, covering reaches into the
interior, so on the same data the two optima are generally different
selections.

### `KCentre()`: near-optimal covering

[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
chooses centres with the deterministic CDSh heuristic ([García-Díaz et
al., 2017](#ref-GarciaDiaz2017); [García-Díaz et al.,
2019](#ref-GarciaDiaz2019)), which typically lands within 1–3.5 % of the
optimum — an order of magnitude tighter than the Gonzalez
2-approximation that
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
gives for this objective.

``` r

centres <- KCentre(eurodist, k = 4L)
labels(eurodist)[centres]
#> [1] "Cologne"    "Copenhagen" "Madrid"     "Rome"
centres
#> 4 centres (6 7 14 19) by CDSh, covering radius <= 1011
```

[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
scores any centre set by its covering radius (lower is better). CDSh
covers at least as tightly as the Gonzalez 2-approximation baseline:

``` r

ff <- FarFirst(eurodist, k = 4L, method = "peripheral")
c(KCentre  = KCentreRadius(eurodist, centres),
  FarFirst = KCentreRadius(eurodist, ff))
#>  KCentre FarFirst 
#>     1011     1209
```

Like [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md),
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
also accepts a `points` coordinate matrix, in which case it never
materialises the full *N × N* matrix — so it can score a selection well
past the size at which the solvers’ distance matrix would fit in memory.

### `ExactKCentre()`: proven covering optimum

For small instances,
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
solves the covering problem to proven optimality with a sequence of
minimum-set-cover integer programs — the covering dual of
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)’s
node-packing programme — warm-started from the
[`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md)
radius and bisected down to the smallest feasible radius. Like
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md),
it uses the **highs** solver.

``` r

res_kc <- ExactKCentre(eurodist, k = 4L, progress = FALSE)
res_kc
#> 4 centres (6 7 14 19) by exact MILP (highs), proven optimal, covering radius = 1011
res_kc$proven      # TRUE  ⟹  radius is the global covering optimum
#> [1] TRUE
```

The covering optimum is sometimes attained by fewer centres (once every
point is covered, extra centres cannot lower the radius); `indices` then
has length below the requested number, and the reported `radius` is
still the proven optimum.
[`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md)
is NP-hard, so — like
[`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md)
— it is a ground-truth reference for small instances, not a scalable
method.

The dispersion and covering optima differ even on this small example:
dispersion selects cities at the rim of the map, while covering pulls
inward to keep every city near a centre.

``` r

disp <- sort(ExactMaxMin(eurodist, k = 4L)$indices)
labels(eurodist)[disp]            # dispersion: pushed to the extremes
#> [1] "Athens"    "Lisbon"    "Milan"     "Stockholm"
labels(eurodist)[sort(res_kc$indices)]  # covering: pulled toward the interior
#> [1] "Cologne"    "Copenhagen" "Madrid"     "Rome"
```

## When to use which method

For **dispersion** (spread the selection; maximise T_(k)):

| Scenario | Recommended |
|----|----|
| Speed matters most, N up to a few thousand | [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) (ensemble default) |
| Deterministic, reproducible refinement | [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) |
| Best quality, [`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility | [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) |
| N \> 46 000 (distance matrix infeasible) | `DropAdd(points = ...)` or `FarFirst(points = ...)` |
| Arbitrary metric with no coordinate embedding | `FarFirst(<column function>, N = ...)` |
| Proven optimum, N ≤ ~ 25–30, **highs** installed | [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) |
| Score a selection’s T_(k) | [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) |

For **covering** (minimise the radius; no point far from a centre):

| Scenario | Recommended |
|----|----|
| Near-optimal covering, fast and deterministic | [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) (CDSh) |
| A quick 2-approximation baseline | `FarFirst(method = "peripheral")` |
| Proven optimum, small N, **highs** installed | [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md) |
| Score a centre set’s covering radius (matrix-free at large N) | `KCentreRadius(points = ...)` |

The dispersion heuristics are complementary:
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) is
O(*N* · *k*) and deterministic — an instant first result.
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) is
deterministic and reproducible (no RNG): ties are broken by smallest
index, so a given instance always yields the same selection, typically
reaching ≈ 99 % of optimal.
[`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) is
stochastic and usually edges out
[`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) on
T_(k); call [`set.seed()`](https://rdrr.io/r/base/Random.html) before it
for a reproducible run.

## Related packages

The CRAN package
[**maximin**](https://CRAN.R-project.org/package=maximin) constructs
continuous space-filling designs by generating new points.

## References

García-Díaz, J., Menchaca-Méndez, R., Menchaca-Méndez, R., Pomares
Hernández, S., Pérez-Sansalvador, J. C., & Lakouari, N. (2019).
Approximation algorithms for the vertex $`k`$-center problem: Survey and
experimental evaluation. *IEEE Access*, *7*, 109228–109245.
<https://doi.org/10.1109/ACCESS.2019.2933875>

García-Díaz, J., Sánchez-Hernández, J., Menchaca-Méndez, R., &
Menchaca-Méndez, R. (2017). When a worse approximation factor gives
better performance: A 3-approximation algorithm for the vertex
$`k`$-center problem. *Journal of Heuristics*, *23*(5), 349–366.
<https://doi.org/10.1007/s10732-017-9345-x>

González, T. F. (1985). Clustering to minimize the maximum intercluster
distance. *Theoretical Computer Science*, *38*, 293–306.
<https://doi.org/10.1016/0304-3975(85)90224-5>

Hochbaum, D. S., & Shmoys, D. B. (1985). A best possible heuristic for
the $`k`$-center problem. *Mathematics of Operations Research*, *10*(2),
180–184. <https://doi.org/10.1287/moor.10.2.180>

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
