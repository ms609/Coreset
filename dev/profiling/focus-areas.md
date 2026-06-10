# Profiling focus areas — MaxMin

Ranked hot paths for the `/profile` rotation. The package sells on speed across
three solvers; the rotation must touch all three on realistically-sized inputs,
each at a **small** and a **large** `k` (= `n`/`m`, the subset size), and on both
the dense-matrix and matrix-free (`points=`) paths where they exist.

| # | Area | Files | Why hot (1 line) | Paths | k regimes | Last profiled | Status |
|---|------|-------|------------------|-------|-----------|---------------|--------|
| 1 | FarFirst greedy + ensemble | `src/maximin.cpp`, `src/maximin_points.cpp`, `R/seed.R`, `R/gonzalez.R` | O(N·n·dim) greedy pass run 3× per default ensemble; redundant MinDist re-scoring | matrix, points, column-oracle | n/N ∈ {0.01 … 0.5} | 2026-06-10 | OPTIMISED |
| 2 | DropAdd tabu search | `src/dropadd.cpp`, `src/dropadd_mf.cpp`, `R/dropadd.R` | `plateau`(5000) iterations, each O(n) record updates; points path O(n·dim) per column | matrix, points | m ∈ {small, ~n/2} | 2026-06-10 | OPTIMISED matrix-seed (11–13×; T-005a); points-path AT-LIMIT (sqrt load-bearing for sum_dist; T-005b) |
| 3 | GraspPR refinement | `src/grasp.cpp`, `R/grasp.R` | extended-improvement local search scans (n−m) candidates × O(m²) sub-matrix; PR over elite pairs | matrix only | m ∈ {small, large} | 2026-06-10 | OPTIMISED (PR/LS 14–16× at large m; T-006) |

Status legend: `NEW`, `PROFILED`, `OPTIMISED`, `AT-LIMIT`, `SKIPPED`.

Notes
- Realistic sizes: FarFirst/DropAdd scale to large n (points path is O(n) memory);
  GraspPR is dense-only and intended for small n (materialises n×n repeatedly).
- "small k" exercises the construction/per-iteration O(n) cost; "large k"
  exercises the O(m²)/recompute terms that dominate at high subset fraction.
- Drivers must stay ≤ ~8 s bare (bound DropAdd/Grasp with `plateau`/`maxIter`).
