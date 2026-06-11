# MaxMin

`MaxMin` provides solvers to select a dispersed subsample of points in
order to maxmize coverage.

The **Max-Min Diversity Problem** (MMDP, the discrete *p*-dispersion
objective) maximises *separation*. It selects $`k`$ elements such that
the minimum distance between any pair of selected elements is as large
as possible; the chosen elements are maximally distinct. This can reward
selections that leave the interior of the set unrepresented.

The **discrete *k*-centre problem** selects $`k`$ elements such that the
largest distance from *any* element to its nearest selected element is
as small as possible. This rewards selections that reach the whole set,
such that each point has a nearby representative; it can pull centres
inward and collapse well-separated modes onto a central compromise.

Both problems can be quickly approximated within a factor of two by the
greedy farthest-first heuristic
([`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md)).

## Solvers

MMDP (max-min dispersion):

| Function | Method | Use |
|----|----|----|
| [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md) | DropAdd tabu search (Porumbel et al. 2011) | ~99%-optimal heuristic |
| [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md) | GRASP with path-relinking metaheuristic (Resende et al. 2010) | Slower but powerful heuristic |
| [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md) | Node-packing integer program (Sayyady & Fathi 2016) | Proven optimum, small `k` (needs `highs`) |

*k*-centre (min-max covering):

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

## Related packages

`MaxMin` selects a *subset of existing elements*. Several established
packages solve neighbouring objectives:

- **k-medoids / k-median** — minimise the *total* (or mean) distance
  from each element to its nearest centre: *average* coverage, or
  representativeness. This is a different objective from any solved
  here, and is well served elsewhere:
  [`cluster::pam()`](https://cran.r-project.org/package=cluster) and
  `clara()` (PAM / FastPAM / FasterPAM),
  [`ClusterR::Cluster_Medoids()`](https://cran.r-project.org/package=ClusterR),
  and [`banditpam`](https://cran.r-project.org/package=banditpam).
  `pam()` holds the full *O(N²)* dissimilarity matrix (and caps at *n* ≤
  65 536); `clara()` samples to scale; `banditpam` is *O(N* log *N)* and
  matrix-free but accepts only coordinate data with built-in metrics (no
  precomputed matrix or custom distance).

- **k-means** ([`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html))
  minimises within-cluster sum of squares around centres that are
  coordinate *means*, not data points, so it is neither a discrete
  k-centre nor a k-medoids solver and applies only to Euclidean
  coordinates.

- [`maximin`](https://cran.r-project.org/package=maximin) constructs
  continuous space-filling designs by generating new points in a
  coordinate box to maximise the minimum inter-point distance.
