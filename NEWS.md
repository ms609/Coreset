# MaxMin 0.0.0.9009 (development)

## Performance

- `Grasp()` is a further 1.3–1.4× faster at `k = 100` across the `plateau`
  ladder, and ~1.1× at small `k`, single-threaded. Selections are unchanged:
  the kernel returns bit-identical indices, objectives, iteration counts and
  path-relinking call counts over the same 1458-cell grid, and still matches
  the pure-R reference exactly.

  Two ideas, both exact. The local search's extended-improvement tie-break
  count is now derived instead of counted: the candidate's new objective
  floor already pins both halves of the pair count (values below a known
  exact floor cannot exist), so the former pair sweep collapses to a couple
  of reads off the summaries the search already maintains. And a
  construction now hands its distance-to-selection array — with witnesses —
  straight to the local search that refines it, which stops re-reading the
  m distance columns the construction just streamed.

# MaxMin 0.0.0.9008 (development)

## Performance

- `Grasp()` runs its search phases on multiple threads when the package is
  built with OpenMP (the default on Linux and Windows). The thread count
  follows the standard `"mc.cores"` option (default 1, single-threaded).
  Random draws happen only on the main thread, in fixed-size batches, and
  batch results merge in a fixed order, so **a given seed returns the
  identical selection at every core count**: `mc.cores` trades wall-clock
  time only, and published results remain reproducible whatever parallelism
  produced them. Locally measured ~5.9× at 8 threads on the canonical
  n = 2000, k = 100 shape, with single-threaded cost unchanged.

- A further ~1.45× single-threaded on the canonical shapes. The local
  search's remaining per-pass pair sweep is replaced by per-member nearest
  summaries maintained across passes with the same witness technique as
  0.0.0.9007 — one structure now supplies the objective floor, its witness
  edge, the extended-improvement pair count and the post-drop rescores —
  and construction folds its candidate scan into the update pass, skipping
  the final update whose results nothing read.

## Breaking changes

- `Grasp()` selections change once relative to 0.0.0.9007: constructions
  draw their random indices from pre-drawn uniform batches through a
  floor-index mapping (bias below 1 part in 2⁵⁰) — the price of the
  core-count-invariant determinism above. Quality at equal search effort and
  at equal time budget is unchanged within noise, measured across seeds,
  shapes and budgets. Seed-reproducibility within 0.0.0.9008 is exact, at
  any core count, and the pure-R reference implementation mirrors the new
  draw protocol bit for bit.

# MaxMin 0.0.0.9007 (development)

## Performance

- `Grasp()` is a further 3.9–4.4× faster at `k = 100` across the whole
  `plateau` ladder (n = 2000 coreset shapes; ~15× below `plateau` 64 and ~9.5×
  at `plateau` 256 relative to 0.0.0.9005), with smaller subset sizes
  unchanged or modestly faster. Selections are unchanged: the kernel returns
  bit-identical indices, objectives, iteration counts and path-relinking call
  counts over the same 1458-cell grid as before, and still matches the pure-R
  reference exactly. Because nothing about the search trajectory moved,
  results computed with earlier versions remain valid.

  Three ideas, all exact. The local search's candidate scan now screens
  through an incrementally maintained nearest-selected-point summary instead
  of rescanning every candidate against the selection: when a candidate's
  stored nearest survives the proposed drop, its post-drop distance is
  already known, and only the few candidates whose nearest is being dropped
  are rescanned. The local search's per-pass bookkeeping folds three
  pair sweeps into one plus a lazily computed exclusion pass. And the
  path-relinking walk carries its candidate lists and nearest summaries
  across steps, recomputing only what a one-element swap can actually
  disturb. At a fixed `maxSeconds` budget the new kernel completes several
  times more GRASP iterations and reaches path relinking where the old one
  ran out of clock; at equal budgets it never returned a worse selection
  across 60 seeded comparisons.

# MaxMin 0.0.0.9006 (development)

- `Grasp()` no longer discards its best-known solution, and is ~2.0 faster.

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
