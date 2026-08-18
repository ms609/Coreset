# Coreset 0.0.0.9002 (development)

- `ExactMaxMin()` gains `nStart`, `graspPlateau` and `dropPlateau`, which set
  how many `Grasp()` restarts seed the search and how long each restart and the
  `DropAdd()` pass run before stopping. The defaults are what the solver
  already used, so no existing call changes. The pool only supplies a lower
  bound, so no setting can alter the optimum reported -- a shallower pool
  simply leaves more thresholds to be searched.

# Coreset 0.0.0.9001 (development)

- `ExactMaxMin()` answers a threshold from the subset it already holds where
  that subset attains it. A feasibility probe returns some `k`-clique rather
  than a best one, and the minimum distance it realises is often well above the
  threshold that was asked for -- on the larger cases, hundreds of attainable
  thresholds above -- so the probes in between were re-finding the same subset.
  The thresholds the search visits are unchanged, and so is every optimum it
  reports. Worth 37% on `tc13_pima` at `k` = 48 and 8% on `tc19_breastcancer`;
  no cell of the benchmark suite was slower, and the sub-second cells sat
  within noise. The gain is a continuous-distance one: where distances are
  integers the witness rarely reaches past the threshold it was asked for.

# Coreset 0.0.0.9000 (development)

First release. `Coreset` selects a representative subset of a fixed candidate
set under an arbitrary distance, solving four discrete location objectives on a
distance matrix, a coordinate matrix, or an on-demand distance-column oracle
(for spaces with no coordinate embedding).

## Max-Min diversity (MMDP / discrete _p_-dispersion)

Maximises the minimum pairwise distance within a subset of size `k`.

- `FarFirst()`: greedy farthest-first selection (Gonzalez 1985), with a choice
  of peripheral seeding strategies, distinct-seed random restarts (`nSeeds`),
  and a robust ensemble default.
- `DropAdd()`: DropAdd tabu search (Porumbel et al. 2011), which can compute
  distances between pairs on the fly rather than needing a complete matrix
  _a priori_.
- `Grasp()`: GRASP with path relinking (Resende et al. 2010), attaining the
  highest `T_k` of the package's heuristics on small to medium instances.
- `ExactMaxMin()`: exact node-packing optimum (Sayyady & Fathi 2016), decided
  by clique search rather than an integer program.

## Max-Mean dispersion

- `MaxMean()`: reinforcement-learning-guided tabu search (Nijimbere et al.
  2020), selecting a subset of unrestricted size that maximises the mean
  pairwise distance.

## Discrete _k_-centre

Minimises the largest distance from any element to its nearest selected centre.

- `KCentre()`: the CDSh covering heuristic (Garcia-Diaz et al. 2017, 2019).
- `ExactKCentre()`: exact minimum-cover optimum (needs `highs`).

## Max-Sum diversity and maximum entropy

- `ExactMaxSum()`: exact solver for the Max-Sum Diversity Problem.
- `MaxEntropy()`: maximum-entropy (maxdet) selection — the mode of a
  determinantal point process — by greedy pivoted-Cholesky selection, and by
  exact enumeration for small instances.

## Scoring and utilities

- `MinDist()`, `MeanDist()` and `KCentreRadius()` score an arbitrary selection
  under the max-min, max-mean and _k_-centre objectives respectively.
- `PickPoint()` exposes the peripheral seed indices directly.
- `DropAdd()` and `Grasp()` accept a `maxCandidates` composable-coreset cap,
  thinning the candidate set with `FarFirst()` before the expensive search and
  mapping the chosen indices back to the original numbering.
- Each solver returns a classed object with `print()`, `format()` and (where
  informative) `summary()` methods giving a terse or detailed report of the
  selection, the achieved objective, and the search effort.
- Solver behaviour is tunable through `options(Coreset.symmetryTolerance = )`,
  which sets how large a rounding discrepancy between `d[i, j]` and `d[j, i]`
  is repaired rather than refused, and `options(Coreset.progress = )`.
