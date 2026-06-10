# MaxMin

`MaxMin` selects a maximally dispersed subset of a fixed set of
candidate items under the Max-Min Diversity Problem (MMDP, the discrete
*p*-dispersion objective): choose $`n`$ items so that the minimum
pairwise distance within the selection is as large as possible.

The package supports precomputed distance matrices (`dist` objects),
matrices of Euclidian coordinates, or lists of elements from which
distances can be calculated.

## Solvers

| Function | Method | Use |
|----|----|----|
| [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) | GRASP with path-relinking metaheuristic (Resende et al. 2010) | Slow but powerful heuristic |
| [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) | Greedy farthest-first (Gonzalez 1985); default best of three random peripheral starts | Fast; matrix, coordinate, or distance-column-oracle input (the last for very large sets with no embedding) |
| [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `n` (needs `highs`) |
| [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md) | Minimum pairwise distance (the objective) | Score a selection |

## Installation

``` r

# install.packages("remotes")
remotes::install_github("ms609/MaxMin")
```

## Related packages

The CRAN package [`maximin`](https://cran.r-project.org/package=maximin)
constructs continuous space-filling designs — it generates new points in
a coordinate box to maximise the minimum inter-point distance.
