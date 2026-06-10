# MaxMin 0.0.0.9002 (development)

## Improvements

- `DropAdd()` now documents that `timeBudgetS` is checked every 256 iterations and
  may overshoot by up to one iteration's worth of computation on large instances.
- `FarFirst()` documents that asymmetric distance matrices are accepted.
- Ensemble functions now attach a `score = NA_real_` attribute on the trivial
  all-points early return (when `m >= N`).
- Integer iteration counters in the C++ kernels changed from `int` to
  `long long` to avoid signed-integer overflow upper bound at extreme `maxIter`
  values.
- Test suite: improved coverage of path relinking (strict improvement), DropAdd
  `secondary` attribute formula, `ExactMaxMin` budget-expiry branch,
  and various weak / vacuous assertions tightened.

# MaxMin 0.0.0.9001 (development)

## Bug fixes

- `Grasp()` no longer crashes at the documented `alpha = 1` (pure greedy).
- `Grasp()` now validates `alpha`, rejecting values outside `[0, 1]`.
- `FarFirst()`, `DropAdd()`, `Grasp()` and `MinDist()` now reject a distance
  matrix containing `NA`/`NaN`/`Inf` instead of silently returning a selection
  with a repeated index; the distance-column oracle path of `FarFirst()`
  likewise rejects a non-self `NA`/`NaN`.
- `MinDist()` now errors on `NA` or duplicate indices in `idx` (previously a
  matrix-path `NA` returned `NA` silently and duplicates scored 0).
- Added defensive guards to the exported C++ kernels reachable via `:::`
  (`MaximinFrom_cpp`, `MaximinFromPoints_cpp`, `DropAdd_points_cpp`).

## Breaking API changes

- `FarFirst()`: the subset-size argument is renamed `n` -> `m`, matching the
  other solvers; and the `seed` argument is renamed `method` (matching
  `MaxMinSeed(method =)`), since it selects a seeding *strategy*, not an RNG
  seed.
- `DropAdd()` and `Grasp()`: the `seed` argument is removed. `DropAdd()`'s was
  a documented no-op (the search is RNG-free); for a reproducible `Grasp()` run,
  call `set.seed()` before `Grasp()`.
- `FarFirst()` now attaches the achieved objective (T_k) as a `score` attribute
  on both bare and ensemble returns, matching `DropAdd()` and `Grasp()`
  (previously this value was discarded on a bare pass).
- Index-order conventions are now documented: `FarFirst()` returns
  farthest-first (greedy) order; `DropAdd()`, `Grasp()` and
  `ExactMaxMin()$indices` return ascending indices. `ExactMaxMin()` continues to
  return a list (the deliberate exception), now documented as such.

# MaxMin 0.0.0.9000 (development)

- Initial release: a tiered toolbox for the Max-Min Diversity Problem (MMDP /
  discrete p-dispersion), extracted from the `FurthestPoint` study package so
  that it can be depended on by CRAN packages (e.g. `TreeSearch`).
- `FarFirst()`: deterministic farthest-first selection from a distance matrix,
  Euclidean coordinates (`points =`), or an on-demand **distance-column oracle**
  (pass a column function as `d`, with `N =`) for spaces with no coordinate
  embedding; the oracle may report the self-distance (length `N`) or omit it
  (length `N - 1`), whichever is simpler to compute. The default `seed` runs a
  best-of-three ensemble of
  `"random_furthest"` starts (whose pivots are drawn with the session RNG; set a
  seed for a reproducible selection, or supply them via the `pivots` argument)
  and keeps the best by `MinDist()`. The deterministic O(*N*) anchors
  (`"centroid"`, `"peripheral"`) and the costlier O(*N*²) anchors (`"diameter"`,
  `"anti_medoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`) are available as opt-in
  strategies.
- `DropAdd()`: DropAdd tabu search (Porumbel et al.
  2011); accepts a `dist` object, distance matrix, or coordinate matrix
  (`points =`).
- `ExactMaxMin()`: exact node-packing optimum (Sayyady & Fathi 2016) via the
  `highs` MILP backend.
- `Grasp()`: GRASP with path relinking (Resende et al. 2010), a
  dense-matrix-only refinement metaheuristic that attains the highest `T_k`
  of the package's methods on small to medium instances.
- `MinDist()`: the k-centre objective (minimum pairwise distance).
- `MaxMinSeed()`: exposes the peripheral seed indices directly.
