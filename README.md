# MaxMin

<!-- badges: start -->
[![R-CMD-check](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/MaxMin/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/MaxMin/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ms609/MaxMin)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`MaxMin` provides solvers to select a dispersed subsample of points in order to
maxmize coverage.

The **Max-Min Diversity Problem** (MMDP, the discrete *p*-dispersion objective)
maximises *separation*. It selects $m$ elements such that the minimum distance
between any pair of selected elements is as large as possible; the chosen
elements are maximally distinct.
This can reward selections that leave the interior of the set unrepresented.

The **discrete *k*-centre problem** selects $k$ elements such that the largest
distance from *any* element to its nearest selected element is as small as
possible. This rewards selections that reach the whole set, such that each point
has a nearby representative; it can pull centres inward and collapse
well-separated modes onto a central compromise.

Both problems can be approximated within a factor of two by the greedy
farthest-first heuristic (`FarFirst()`).


## Solvers

MMDP (max-min dispersion):

| Function | Method | Use |
|---|---|---|
| `DropAdd()` | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| `Grasp()` |  GRASP with path-relinking metaheuristic (Resende et al. 2010) | Slow but powerful heuristic |
| `FarFirst()` | Greedy farthest-first (Gonzalez 1985); default best of three random peripheral starts | Fast; matrix, coordinate, or distance-column-oracle input (the last for very large sets with no embedding) |
| `ExactMaxMin()` | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `n` (needs `highs`) |

*k*-centre (min-max covering):

| Function | Method | Use |
|---|---|---|
| `KCentre()` | CDSh heuristic (García-Díaz et al. 2017, 2019) | ~1–3.5% of optimum at $O(N^2 \log N)$, far tighter than the Gonzalez 2-approximation; never worse than it |
| `ExactKCentre()` | Min-cover integer program (the covering dual of `ExactMaxMin()`) | Proven optimum, small `k` (needs `highs`) |

Solvers support precomputed distance matrices (`dist` objects),
matrices of Euclidian coordinates, or lists of elements from which distances
can be calculated.

`MinDist()` returns the minimum pairwise distance within a selection (the MMDP
objective) and `KCentreRadius()` returns its covering radius (the *k*-centre 
objective).


## Installation

```r
# install.packages("remotes")
remotes::install_github("ms609/MaxMin")
```

## Related packages

The CRAN package [`maximin`](https://cran.r-project.org/package=maximin)
constructs continuous space-filling designs — it generates new points in a
coordinate box to maximise the minimum inter-point distance.

