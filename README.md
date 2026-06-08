# MaxMin

<!-- badges: start -->
[![R-CMD-check](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/MaxMin/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/MaxMin)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`MaxMin` selects a maximally dispersed subset of a fixed set of candidate items
under the **Max-Min Diversity Problem** (MMDP, the discrete *p*-dispersion
objective): choose `n` items so that the *minimum* pairwise distance within the
selection is as large as possible.

It is a small, dependency-light toolbox (`Imports: Rcpp, stats`) that operates
on three kinds of input:

- a **distance matrix** (or `dist` object);
- **Euclidean coordinates**, without ever materialising the distance matrix; or
- an **on-demand distance-column oracle**, for spaces with no coordinate
  embedding (e.g. phylogenetic trees), where the distance matrix would be too
  large to hold.

## Solvers

| Function | Method | Use |
|---|---|---|
| `Gonzalez()` | Greedy farthest-first (Gonzalez 1985), with peripheral seeding strategies and an `"ensemble"` default | Fast; matrix, coordinate, or distance-column-oracle input (the last for very large sets with no embedding) |
| `DropAddTS()` / `DropAddTSPoints()` | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| `ExactMaxMin()` | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `n` (needs `highs`) |
| `PolishSelection()` | Critical-edge 1-swap local search | Refine any selection |
| `TkScore()` | Minimum pairwise distance (the objective) | Score a selection |

## Not to be confused with `maximin`

The CRAN package [`maximin`](https://cran.r-project.org/package=maximin) (Sun &
Gramacy) constructs continuous **space-filling designs** — it *generates new
points* in a coordinate box to maximise the minimum inter-point distance.
`MaxMin` instead *selects a subset* from a *fixed* candidate set under an
arbitrary distance, a combinatorial problem on a different footing.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ms609/MaxMin")
```
