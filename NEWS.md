# MaxMin 0.0.0.9000 (development)

- Initial release: a tiered toolbox for the Max-Min Diversity Problem (MMDP /
  discrete p-dispersion), extracted from the `FurthestPoint` study package so
  that it can be depended on by CRAN packages (e.g. `TreeSearch`).
- `Gonzalez()`: deterministic farthest-first selection from a distance matrix,
  Euclidean coordinates (`points =`), or an on-demand **distance-column oracle**
  (pass a column function as `d`, with `N =`) for spaces with no coordinate
  embedding. The default `seed` runs a best-of-five ensemble of cheap O(*N*)
  seeds — `"centroid"`, `"peripheral"`, and three reproducible
  `"random_furthest"` starts whose pivots are drawn from a fixed internal seed
  (override or extend them with the `pivots` argument) — and keeps the best by
  `TkScore()`. Costlier O(*N*²) anchors (`"diameter"`, `"anti_medoid"`,
  `"medoid"`, `"rowsum"`, `"rownorm"`) are available as opt-in strategies.
- `DropAddTS()` / `DropAddTSPoints()`: DropAdd tabu search (Porumbel et al.
  2011), matrix and matrix-free coordinate paths.
- `ExactMaxMin()`: exact node-packing optimum (Sayyady & Fathi 2016) via the
  `highs` MILP backend.
- `GraspPR()`: GRASP with path relinking (Resende et al. 2010), a
  dense-matrix-only refinement metaheuristic that attains the highest `T_k`
  of the package's methods on small to medium instances.
- `PolishSelection()`: critical-edge-anchored 1-swap local search.
- `TkScore()`: the k-centre objective (minimum pairwise distance).
- `MaxMinSeed()`: exposes the peripheral seed indices directly.
