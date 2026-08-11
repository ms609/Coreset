# MaxMin 0.0.0.9005 (development)

- `DropAdd()` gains a **distance-column oracle** path, parallel to
  `FarFirst()`'s: pass a function as `d` (with `N`) and it is called as `d(i)`
  to obtain the distances from element `i` to the others, one column at a time.
  This suits metrics with neither a stored matrix nor a coordinate embedding --
  on-demand tree-to-tree distances, say -- and never materialises the `N x N`
  matrix, so memory is `O(N)`. The search itself is unchanged: on a symmetric
  oracle the whole drop-add trajectory, and hence the selection, `score` and
  `secondary`, are bit-identical to the matrix path given the same `seed`.
  Two things differ, both documented under `?DropAdd`: the max-row-sum warm
  start would cost all `N` columns, so an `O(N)` two-sweep peripheral seed is
  substituted (override with `seed=`); and `maxCandidates` thinning is
  unavailable, warning rather than silently thinning when the cap binds.
  Cost is counted in oracle calls (`k` to construct, two per iteration), so
  `plateau` and `maxSeconds` set the distance bill directly -- see
  `?DropAdd`'s example for the closure shape that keeps it cheap.

# MaxMin 0.0.0.9004 (development)

- New `MaxEntropy()`: maximum-entropy (maxdet) subset selection -- choose the
  `k`-subset maximising `log det K_S` for a radial-basis similarity kernel built
  from the distances (Shewry & Wynn 1987; the mode of a determinantal point
  process, Kulesza & Taskar 2012). The kernel is repaired to positive-
  semidefinite form (eigen-clip by default); selection is by greedy pivoted
  Cholesky, with exact enumeration where `choose(n, k)` is small. A redundant
  point adds zero volume and is never co-selected, so the objective is exactly
  density-blind. Returns a `"MaxEntropySelection"` carrying the retained
  log-determinant and the repair magnitude, with its own `format()`/`print()`.

- New `ExactMaxSum()`: exact solver for the Max-Sum Diversity Problem (the
  "maximum diversity problem" -- select the `k`-subset of greatest *total*
  pairwise distance), the max-sum counterpart of `ExactMaxMin()`. Uses the
  Kuo--Glover--Dhir per-node MILP linearisation (`highs`) with a multi-start
  1-swap local-search incumbent as floor and time-limited fallback. Returns a
  `"MaxSumSelection"` with its own `format()`/`print()` reporting the achieved
  total distance.

- `DropAdd()` and `Grasp()` gain a `maxCandidates` argument: a composable
  coreset (Indyk et al. 2014; Aghamolaei et al. 2015) that thins the candidates
  to a farthest-first subset before the heavy search, letting the solvers run at
  scales where they were previously intractable. The coreset uses the
  deterministic `"peripheral"` seed (no session RNG is drawn) and the chosen
  indices are mapped back to the original numbering. Thinning is **on by
  default** (`DropAdd()`: `46340`; `Grasp()`: `2000`) and emits a warning when it
  binds; pass `maxCandidates = 0` to run on the full problem as before. A request
  for more points than the cap (`k > maxCandidates`) is now an error.

- `FarFirst()`'s default `nSeeds` is reduced from `8` to `3`. The random-furthest
  restart gain curve bends early (a knee at roughly three to four starts across
  benchmarks), so three starts capture most of the benefit at well under half the
  cost. This changes the default selection for a given `set.seed()`; pass
  `nSeeds = 8` for the previous default, or `DropAdd()` for higher-quality results.

- `DropAdd()` gains a `seed` argument: an optional 1-based start index that
  overrides the construction's default warm-start (the max-row-sum point on the
  `d` path, the centroid-peripheral point on the `points` path), mirroring
  `FarFirst()`'s integer `strategy`. `NULL` (the default) keeps the method's own
  seed. It is not supported when candidate thinning binds (pass
  `maxCandidates = 0`).

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
