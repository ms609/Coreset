# Changelog

## Coreset 1.0.0.9004 (development)

- [`DropAdd()`](https://ms609.github.io/Coreset/dev/reference/DropAdd.md)
  on a distance matrix no longer crashes beyond 46,340 points, where its
  column offsets overflowed an `int`;
  [`ExactMaxMin()`](https://ms609.github.io/Coreset/dev/reference/ExactMaxMin.md)’s
  warm start calls it, so this also affected the exact solver at that
  size.

- [`ExactMaxMin()`](https://ms609.github.io/Coreset/dev/reference/ExactMaxMin.md)
  reports an `upper` bound on the optimum, and can spend `boundSeconds`
  bracketing the optimum from above before its search, so an unproven
  result still bounds its distance from the optimum. Passing that bound
  back as `upper`, with the selection as `warmStart`, resumes the
  search.

- [`PickPoint()`](https://ms609.github.io/Coreset/dev/reference/PickPoint.md)
  supports the `"peripheral"` and `"random_furthest"` strategies,
  including where `d` is a function, returning `nSeeds` seeds.

- [`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md)
  supports `strategy = "random_furthest"` when `d` is a function.

- [`MaxEntropy()`](https://ms609.github.io/Coreset/dev/reference/MaxEntropy.md)
  is roughly 3-12x faster at large `n`.

## Coreset 1.0.0 (2026-09-09)

CRAN release: 2026-09-17

`Coreset` selects a representative subset of a fixed candidate set under
an arbitrary distance, solving four discrete location objectives on a
distance matrix, a coordinate matrix, or an on-demand distance-column
oracle.

### Max-Min diversity (MMDP / discrete *p*-dispersion)

Maximises the minimum pairwise distance within a subset of size `k`.

- [`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md):
  greedy farthest-first selection (Gonzalez 1985), with a choice of
  peripheral seeding strategies, distinct-seed random restarts
  (`nSeeds`), and a robust ensemble default.
- [`DropAdd()`](https://ms609.github.io/Coreset/dev/reference/DropAdd.md):
  DropAdd tabu search (Porumbel et al. 2011), which can compute
  distances between pairs on the fly rather than needing a complete
  matrix *a priori*.
- [`Grasp()`](https://ms609.github.io/Coreset/dev/reference/Grasp.md):
  GRASP with path relinking (Resende et al. 2010), attaining the highest
  `T_k` of the package’s heuristics on small to medium instances.
- [`ExactMaxMin()`](https://ms609.github.io/Coreset/dev/reference/ExactMaxMin.md):
  exact node-packing optimum (Sayyady & Fathi 2016), decided by clique
  search.

### Max-Mean dispersion

- [`MaxMean()`](https://ms609.github.io/Coreset/dev/reference/MaxMean.md):
  reinforcement-learning-guided tabu search (Nijimbere et al. 2020),
  selecting a subset of unrestricted size that maximises the mean
  pairwise distance.

### Discrete *k*-centre

Minimises the largest distance from any element to its nearest selected
centre.

- [`KCentre()`](https://ms609.github.io/Coreset/dev/reference/KCentre.md):
  the CDSh covering heuristic (Garcia-Diaz et al. 2017, 2019).
- [`ExactKCentre()`](https://ms609.github.io/Coreset/dev/reference/ExactKCentre.md):
  exact minimum-cover optimum.

### Max-Sum diversity and maximum entropy

- [`ExactMaxSum()`](https://ms609.github.io/Coreset/dev/reference/ExactMaxSum.md):
  exact solver for the Max-Sum Diversity Problem (requires ‘highs’).
- [`MaxEntropy()`](https://ms609.github.io/Coreset/dev/reference/MaxEntropy.md):
  maximum-entropy (maxdet) selection — the mode of a determinantal point
  process — by greedy pivoted-Cholesky selection, and by exact
  enumeration for small instances.

### Scoring and utilities

- [`MinDist()`](https://ms609.github.io/Coreset/dev/reference/MinDist.md),
  [`MeanDist()`](https://ms609.github.io/Coreset/dev/reference/MeanDist.md)
  and
  [`KCentreRadius()`](https://ms609.github.io/Coreset/dev/reference/KCentreRadius.md)
  score an arbitrary selection under the max-min, max-mean and
  *k*-centre objectives respectively.
- [`PickPoint()`](https://ms609.github.io/Coreset/dev/reference/PickPoint.md)
  exposes the peripheral seed indices directly.
- [`DropAdd()`](https://ms609.github.io/Coreset/dev/reference/DropAdd.md)
  and
  [`Grasp()`](https://ms609.github.io/Coreset/dev/reference/Grasp.md)
  accept a `maxCandidates` composable-coreset cap, thinning the
  candidate set with
  [`FarFirst()`](https://ms609.github.io/Coreset/dev/reference/FarFirst.md)
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
