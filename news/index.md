# Changelog

## Coreset 0.0.0.9000 (development)

First release. `Coreset` selects a representative subset of a fixed
candidate set under an arbitrary distance, solving four discrete
location objectives on a distance matrix, a coordinate matrix, or an
on-demand distance-column oracle (for spaces with no coordinate
embedding).

### Max-Min diversity (MMDP / discrete *p*-dispersion)

Maximises the minimum pairwise distance within a subset of size `k`.

- [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md):
  greedy farthest-first selection (Gonzalez 1985), with a choice of
  peripheral seeding strategies, distinct-seed random restarts
  (`nSeeds`), and a robust ensemble default.
- [`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md):
  DropAdd tabu search (Porumbel et al. 2011), which can compute
  distances between pairs on the fly rather than needing a complete
  matrix *a priori*.
- [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md): GRASP
  with path relinking (Resende et al. 2010), attaining the highest `T_k`
  of the package’s heuristics on small to medium instances.
- [`ExactMaxMin()`](https://ms609.github.io/Coreset/reference/ExactMaxMin.md):
  exact node-packing optimum (Sayyady & Fathi 2016), decided by clique
  search rather than an integer program.

### Max-Mean dispersion

- [`MaxMean()`](https://ms609.github.io/Coreset/reference/MaxMean.md):
  reinforcement-learning-guided tabu search (Nijimbere et al. 2020),
  selecting a subset of unrestricted size that maximises the mean
  pairwise distance.

### Discrete *k*-centre

Minimises the largest distance from any element to its nearest selected
centre.

- [`KCentre()`](https://ms609.github.io/Coreset/reference/KCentre.md):
  the CDSh covering heuristic (Garcia-Diaz et al. 2017, 2019).
- [`ExactKCentre()`](https://ms609.github.io/Coreset/reference/ExactKCentre.md):
  exact minimum-cover optimum (needs `highs`).

### Max-Sum diversity and maximum entropy

- [`ExactMaxSum()`](https://ms609.github.io/Coreset/reference/ExactMaxSum.md):
  exact solver for the Max-Sum Diversity Problem.
- [`MaxEntropy()`](https://ms609.github.io/Coreset/reference/MaxEntropy.md):
  maximum-entropy (maxdet) selection — the mode of a determinantal point
  process — by greedy pivoted-Cholesky selection, and by exact
  enumeration for small instances.

### Scoring and utilities

- [`MinDist()`](https://ms609.github.io/Coreset/reference/MinDist.md),
  [`MeanDist()`](https://ms609.github.io/Coreset/reference/MeanDist.md)
  and
  [`KCentreRadius()`](https://ms609.github.io/Coreset/reference/KCentreRadius.md)
  score an arbitrary selection under the max-min, max-mean and
  *k*-centre objectives respectively.
- [`PickPoint()`](https://ms609.github.io/Coreset/reference/PickPoint.md)
  exposes the peripheral seed indices directly.
- [`DropAdd()`](https://ms609.github.io/Coreset/reference/DropAdd.md)
  and [`Grasp()`](https://ms609.github.io/Coreset/reference/Grasp.md)
  accept a `maxCandidates` composable-coreset cap, thinning the
  candidate set with
  [`FarFirst()`](https://ms609.github.io/Coreset/reference/FarFirst.md)
  before the expensive search and mapping the chosen indices back to the
  original numbering.
- Each solver returns a classed object with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`format()`](https://rdrr.io/r/base/format.html) and (where
  informative) [`summary()`](https://rdrr.io/r/base/summary.html)
  methods giving a terse or detailed report of the selection, the
  achieved objective, and the search effort.
- Solver behaviour is tunable through
  `options(Coreset.symmetryTolerance = )`, which sets how large a
  rounding discrepancy between `d[i, j]` and `d[j, i]` is repaired
  rather than refused, and `options(Coreset.progress = )`.
