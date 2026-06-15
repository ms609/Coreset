# MaxMin

`MaxMin` implements algorithms that select a dispersed subsample of
points that maximizes coverage of a larger set, under one of two
objectives:

The **Max-Min Diversity Problem** (MMDP, the discrete *p*-dispersion
objective) selects $`k`$ elements such that the minimum distance between
any pair of selected elements is as large as possible; the chosen
elements are maximally separated.

The **discrete *k*-centre problem** selects $`k`$ elements such that the
maximum distance from any element in the original set to a selected
element is as small as possible.

## Solvers

The greedy farthest-first heuristic
[`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)
serves as a quick approximate solver to both problems, providing a
solution that is guaranteed to be within a factor of two of the optimum
(and typically much closer).

Better approximations can be accomplished, at the cost of longer runs,
by specialized solvers.

### MMDP (max-min dispersion)

| Function | Method | Use |
|----|----|----|
| [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) | GRASP with path-relinking metaheuristic (Resende et al. 2010) | Slower but powerful heuristic |
| [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `k` (needs `highs`) |

### *k*-centre (min-max covering)

| Function | Method | Use |
|----|----|----|
| [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md) | CDSh heuristic (García-Díaz et al. 2017, 2019) | ~1–3.5% of optimum at $`O(N^2 \log N)`$, typically far tighter than [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md) |
| [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md) | Min-cover integer program | Proven optimum, small `k` (needs `highs`) |

Solvers support precomputed distance matrices (`dist` objects), matrices
of Euclidian coordinates, or lists of elements from which distances can
be calculated.

[`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md)
returns the minimum pairwise distance within a selection (the MMDP
objective) and
[`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md)
returns its covering radius (the *k*-centre objective).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("ms609/MaxMin")
```

## Related problems

`MaxMin` selects a subset from a given set of elements. Several
established packages solve neighbouring objectives:

- **k-medoids / k-median** selects elements that minimize the mean
  distance from each element to its nearest centre. Implementations
  include:

  - [`cluster::pam()`](https://cran.r-project.org/package=cluster):
    generates the full *O(N²)* dissimilarity matrix, and hence caps at
    *n* ≤ 65 536;
  - [`banditpam`](https://cran.r-project.org/package=banditpam), a
    matrix-free $`O(N \log N)`$ implementation restricted to coordinate
    data;
  - [`cluster::clara()`](https://rdrr.io/pkg/cluster/man/clara.html)
    (PAM / FastPAM / FasterPAM);
  - [`ClusterR::Cluster_Medoids()`](https://cran.r-project.org/package=ClusterR).

- **k-means** ([`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html))
  selects elements so as to minimize the within-cluster sum of squares
  around centres that are coordinate means, not data points; as such, it
  applies only to Euclidean coordinates. k-means++
  ([`TreeDist::KMeansPP()`](https://ms609.github.io/TreeDist/reference/KMeansPP.html))
  initializes its selection using D²-weighted seeding, a randomized
  relative of
  [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)’s
  farthest-first traversal.

- [`maximin`](https://cran.r-project.org/package=maximin) solves the
  related design problem of adding *new* points at positions that
  maximize the minimum inter-point distance.
