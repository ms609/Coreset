# Coreset

<!-- badges: start -->
[![R-CMD-check](https://github.com/ms609/Coreset/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/Coreset/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/Coreset/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/Coreset)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`Coreset` implements algorithms for discrete diversity, dispersion, and coverage
subset selection: choosing a representative subsample of points under one of
several objectives.

## Solvers

Solvers support precomputed distance matrices (`dist` objects),
matrices of Euclidian coordinates, or lists of elements from which distances
can be calculated. Each solver is accompanied by a function that evaluates
its objective against an arbitrary selection of elements.

### MMDP (max-min dispersion)

The Max-Min Diversity Problem (MMDP, the discrete *p*-dispersion objective)
selects $k$ elements such that the minimum distance between any pair of selected elements is as large as possible; the chosen elements are maximally separated.

| Function | Method | Use |
|---|---|---|
| `DropAdd()` | DropAdd tabu search | ~99%-optimal heuristic |
| `Grasp()` |  GRASP with path-relinking metaheuristic | Slower but powerful heuristic |
| `ExactMaxMin()` | Node-packing integer program | Proven optimum, small `k` |

### Max-sum dispersion (maximum diversity)

The Max-Sum dispersion problem selects $k$ elements such that the summed
distance between all pairs of selected elements is as large as possible.

| Function | Method | Use |
|---|---|---|
| `ExactMaxSum()` | Per-node MILP linearisation, floored by multi-start local search | Proven optimum for total pairwise distance, small `k` (needs `highs`) |

### Max-mean dispersion

The analogous Max-Mean dispersion problem chooses a number of elements
so as to maximize the mean distance between selected pairs.

| Function | Method | Use |
|---|---|---|
| `MaxMean()` | Reinforcement-learning tabu search | Maximize mean pairwise distance; subset size free; signed distances supported |

### *k*-centre (min-max covering)

The discrete *k*-centre problem selects $k$ elements such that the maximum
distance from any element in the original set to a selected element is as small
as possible.

| Function | Method | Use |
|---|---|---|
| `KCentre()` | Critical Dominating Set heuristic | ~1–3.5% of optimum at $O(N^2 \log N)$, typically far tighter than `FarFirst()` |
| `ExactKCentre()` | Min-cover integer program | Proven optimum, small `k` (needs `highs`) |

### Maximum-entropy (maxdet) selection

The maximum entropy problem selects the $k$ elements that contain the highest
amount of information about the original set: a minimally redundant pick.

| Function | Method | Use |
|---|---|---|
| `MaxEntropy()` | Greedy pivoted-Cholesky selection, with exact enumeration for small instances | Maximize the log-determinant (spanned volume) of a similarity kernel; density-blind |


## Installation

```r
# install.packages("remotes")
remotes::install_github("ms609/Coreset")
```

## Related problems

`Coreset` selects a subset from a given set of elements.
Several established packages solve neighbouring objectives:

- **k-medoids / k-median** selects elements that minimize the mean distance from
  each element to its nearest centre.
  Implementations include:
  * [`cluster::pam()`](https://cran.r-project.org/package=cluster): generates
  the full *O(N²)* dissimilarity matrix, and hence caps at *n* ≤ 65 536;
  * [`banditpam`](https://cran.r-project.org/package=banditpam), a matrix-free
  $O(N \log N)$ implementation restricted to coordinate data;
  * `cluster::clara()`  (PAM / FastPAM / FasterPAM);
  * [`ClusterR::Cluster_Medoids()`](https://cran.r-project.org/package=ClusterR).

- **k-means** ([`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)) selects
  elements so as to minimize the within-cluster sum of squares around centres
  that are coordinate means, not data points; as such, it applies only to
  Euclidean coordinates. k-means++ 
  ([`TreeDist::KMeansPP()`](
   https://ms609.github.io/TreeDist/reference/KMeansPP.html)) initializes its
   selection using D²-weighted seeding, a randomized relative of `FarFirst()`'s
   farthest-first traversal.

- [`maximin`](https://cran.r-project.org/package=maximin) solves the related
  design problem of adding *new* points at positions that maximize the minimum
  inter-point distance.
