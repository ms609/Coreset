# MaxMin 0.0.0.9008 (development)

- `MaxMean()` solves the Max-Mean Dispersion Problem.
- `MeanDist()` scores an arbitrary selection under the max-mean objective.

# MaxMin 0.0.0.9007 (development)

- `FarFirst()` parallelized and optimized.

# MaxMin 0.0.0.9006 (development)

- `Grasp()` no longer discards its best-known solution, is ~10.0× faster,
  and may be parallelized.

# MaxMin 0.0.0.9005 (development)

- `DropAdd()` can now compute distances between pairs on the fly, rather than
  needing a complete matrix _a priori_.
  
# MaxMin 0.0.0.9004 (development)

- New `MaxEntropy()`: maximum-entropy (maxdet) subset selection.

- New `ExactMaxSum()`: exact solver for the Max-Sum Diversity Problem.

- `DropAdd()` and `Grasp()` gain a `maxCandidates` argument.

- `FarFirst()`'s default `nSeeds` is reduced from `8` to `3`.

- `DropAdd()` gains a `seed` argument.

# MaxMin 0.0.0.9003 (development)

- `KCentre()` / `ExactKCentre()` solve the k-centre problem.

# MaxMin 0.0.0.9002 (development)

## Improvements

- The solver results now print legibly. `FarFirst()`, `DropAdd()` and `Grasp()`
  return a `"MaxMinSelection"` object and `ExactMaxMin()` a `"MaxMinExact"`
  object, each with a `print()`/`format()` method giving a one-line summary
  (size, selected indices, algorithm and -- for a `FarFirst()` ensemble -- the
  winning strategy, plus the achieved `T_k`). A `summary()` method adds a
  multi-line report: the per-strategy `T_k` table for a `FarFirst()` ensemble,
  the secondary objective and search effort for `DropAdd()`/`Grasp()`, and the
  instance, objective and proof status for `ExactMaxMin()`. The objects are
  otherwise unchanged: a `MaxMinSelection` is still the integer index vector (it
  indexes a matrix or coordinate set directly), and a `MaxMinExact` is still the
  list with `$indices`, `$objective`, ....

- `FarFirst()` gains an `nSeeds` argument: a distinct-seed random restart that
  draws random pivots, collects each one's furthest-point seed de-duplicated
  until `nSeeds` *distinct* seeds are found (or the reachable pool is exhausted),
  runs Gonzalez from each and returns the best `T_k`. It is the "give a count,
  not a list" counterpart to `pivots` -- never wasting a Gonzalez pass on a
  duplicate seed -- and overrides `method`/`pivots` when supplied. Reproducible
  under `set.seed()`; unsupported on the distance-column oracle path.

- `ExactMaxMin()` is substantially faster and scales to larger instances.
- `DropAdd()` now documents that `maxSeconds` is checked every 256 iterations and
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
  `PickPoint(method =)`), since it selects a seeding *strategy*, not an RNG
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
  discrete p-dispersion).
- `FarFirst()`: deterministic farthest-first selection from a distance matrix,
  Euclidean coordinates (`points =`), or an on-demand **distance-column oracle**
  (pass a column function as `d`, with `N =`) for spaces with no coordinate
  embedding; the oracle may report the self-distance (length `N`) or omit it
  (length `N - 1`), whichever is simpler to compute. The default `seed` runs a
  best-of-three ensemble of
  `"random_furthest"` starts (whose pivots are drawn with the session RNG; set a
  seed for a reproducible selection, or supply them via the `pivots` argument)
  and keeps the best by `MinDist()`. The deterministic O(*N*) anchors
  (`"anti_centroid"`, `"peripheral"`) and the costlier O(*N*²) anchors (`"diameter"`,
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
- `PickPoint()`: exposes the peripheral seed indices directly.
