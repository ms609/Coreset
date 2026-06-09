# MaxMin

`MaxMin` selects a maximally dispersed subset of a fixed set of
candidate items under the Max-Min Diversity Problem (MMDP, the discrete
*p*-dispersion objective): choose $`n`$ items so that the minimum
pairwise distance within the selection is as large as possible.

This dependency-light toolbox operates on three types of input:

- a distance matrix (or `dist` object);
- Euclidean coordinates, without ever materialising the distance matrix;
  or
- an on-demand distance-column oracle, for spaces with no coordinate
  embedding (e.g. phylogenetic trees), where the distance matrix would
  be too large to store in memory.

## Solvers

| Function | Method | Use |
|----|----|----|
| [`Gonzalez()`](https://ms609.github.io/MaxMin/reference/Gonzalez.md) | Greedy farthest-first (Gonzalez 1985); default best-of-three ensemble of reproducible random-furthest starts (deterministic anchors such as centroid/peripheral are opt-in) | Fast; matrix, coordinate, or distance-column-oracle input (the last for very large sets with no embedding) |
| [`DropAddTS()`](https://ms609.github.io/MaxMin/reference/DropAddTS.md) / [`DropAddTSPoints()`](https://ms609.github.io/MaxMin/reference/DropAddTSPoints.md) | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
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
