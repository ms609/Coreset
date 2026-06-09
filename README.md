# MaxMin

<!-- badges: start -->
[![R-CMD-check](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/MaxMin/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/MaxMin)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`MaxMin` selects a maximally dispersed subset of a fixed set of candidate items
under the Max-Min Diversity Problem (MMDP, the discrete *p*-dispersion
objective): choose $n$ items so that the minimum pairwise distance within the
selection is as large as possible.

The package supports precomputed distance matrices (or `dist` objects), or 
matrices of Euclidian coordinates, or lists of elements from which distances
can be calculated.

## Solvers

| Function | Method | Use |
|---|---|---|
| `Gonzalez()` | Greedy farthest-first (Gonzalez 1985); default best-of-three ensemble of reproducible random-furthest starts (deterministic anchors such as centroid/peripheral are opt-in) | Fast; matrix, coordinate, or distance-column-oracle input (the last for very large sets with no embedding) |
| `DropAdd()` / `DropAddPoints()` | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| `ExactMaxMin()` | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `n` (needs `highs`) |
| `MinDist()` | Minimum pairwise distance (the objective) | Score a selection |


## Installation

```r
# install.packages("remotes")
remotes::install_github("ms609/MaxMin")
```

## Related packages

The CRAN package [`maximin`](https://cran.r-project.org/package=maximin)
constructs continuous space-filling designs — it generates new points in a
coordinate box to maximise the minimum inter-point distance.
