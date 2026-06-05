# MaxMin 0.0.0.9000 (development)

- Initial release: a tiered toolbox for the Max-Min Diversity Problem (MMDP /
  discrete p-dispersion), extracted from the `FurthestPoint` study package so
  that it can be depended on by CRAN packages (e.g. `TreeSearch`).
- `Gonzalez()`: deterministic farthest-first selection on a distance matrix or
  Euclidean coordinates, with a `seed` argument choosing a peripheral seeding
  strategy (`"diameter"`, `"antimedoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`)
  and a robust `"ensemble"` default that keeps the best by `TkScore()`.
- `GonzalezColumn()`: matrix-free farthest-first from an on-demand
  distance-column oracle, for spaces with no coordinate embedding.
- `DropAddTS()` / `DropAddTSPoints()`: DropAdd tabu search (Porumbel et al.
  2011), matrix and matrix-free coordinate paths.
- `ExactMaxMin()`: exact node-packing optimum (Sayyady & Fathi 2016) via the
  `highs` MILP backend.
- `PolishSelection()`: critical-edge-anchored 1-swap local search.
- `TkScore()`: the k-centre objective (minimum pairwise distance).
- `MaxMinSeed()`: exposes the peripheral seed indices directly.
