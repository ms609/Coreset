# MaxMin: Maximum-Minimum Diversity, Discrete k-Centre, and Max-Mean Dispersion Subset Selection

Selects a representative subset of a fixed candidate set.

## Max-Min diversity solvers

The Max-Min Diversity Problem (MMDP) maximises the minimum pairwise
distance within a subset (the discrete *p*-dispersion objective).

- [`FarFirst()`](https://ms609.github.io/MaxMin/reference/FarFirst.md):

  Greedy farthest-first selection from a distance matrix, a coordinate
  matrix, or a distance-column oracle (for spaces with no coordinate
  embedding), with a choice of peripheral seeding strategies and a
  robust ensemble default.

- [`DropAdd()`](https://ms609.github.io/MaxMin/reference/DropAdd.md):

  DropAdd tabu search heuristic.

- [`Grasp()`](https://ms609.github.io/MaxMin/reference/Grasp.md):

  GRASP with path relinking.

- [`ExactMaxMin()`](https://ms609.github.io/MaxMin/reference/ExactMaxMin.md):

  Exact node-packing optimum (needs highs).

## Max-Mean dispersion solver

The Max-Mean Dispersion Problem selects a subset of a size that
maximises the mean pairwise distance.

- [`MaxMean()`](https://ms609.github.io/MaxMin/reference/MaxMean.md):

  Reinforcement-learning tabu search.

## k-centre solvers

The discrete *k*-centre problem minimises the largest distance from any
element to its nearest selected element ('centre').

- [`KCentre()`](https://ms609.github.io/MaxMin/reference/KCentre.md):

  CDSh covering heuristic.

- [`ExactKCentre()`](https://ms609.github.io/MaxMin/reference/ExactKCentre.md):

  Exact minimum-cover optimum (needs highs).

## Scoring

- [`MinDist()`](https://ms609.github.io/MaxMin/reference/MinDist.md):

  Minimum pairwise distance (the max-min objective).

- [`MeanDist()`](https://ms609.github.io/MaxMin/reference/MeanDist.md):

  Mean pairwise dispersion (the max-mean objective).

- [`KCentreRadius()`](https://ms609.github.io/MaxMin/reference/KCentreRadius.md):

  Covering radius (the k-centre objective).

## Relation to maximin

Not to be confused with the CRAN package maximin, which constructs
continuous space-filling designs by generating *new* points in a
coordinate region to maximise the minimum inter-point distance.

## See also

Useful links:

- <https://ms609.github.io/MaxMin/>

- Report bugs at <https://github.com/ms609/MaxMin/issues>

## Author

**Maintainer**: Martin R. Smith <martin.smith@durham.ac.uk>
([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder\]

Authors:

- Martin R. Smith <martin.smith@durham.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder\]
