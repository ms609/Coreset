# Agent notes for `Coreset`

`Coreset` is a small, dependency-light (`Imports: cli, Rcpp, Rdpack, stats`) solver library
for discrete diversity, dispersion, and coverage subset-selection problems (Max-Min /
Max-Sum / Max-Mean diversity, discrete k-centre, and max-entropy selection).
**Keep it CRAN-clean and dependency-light** — `highs` (exact solver) stays in
 `Suggests` behind a `requireNamespace()` guard.

## Scope

The MMDP algorithms — Gonzalez farthest-first (matrix / coordinate / 
column-oracle paths) and its seeding strategies, DropAdd tabu search, 
the exact node-packing solver, the polish local search, and the `MinDist`
objective; and the k-center algorithms.

## Conventions (ms609 house style)

- Object naming: exported functions and S3-free helpers in `CamelCase`; internal
  helpers prefixed `.` (e.g. `.MaximinFrom`). The `.lintr` enforces
  CamelCase/camelCase.
- Never use super-assignment (`<<-`). For mutable local caches use
  `new.env(parent = emptyenv())` + the `%||%` helper.
- After any user-visible change, update `NEWS.md` and bump the `.900X` dev
  suffix in `DESCRIPTION`.
- All new/changed code needs test coverage (codecov gates PRs). Cover happy
  paths, error branches, and edge cases (`n == 0`, `n == 1`, `n >= N`).
- Benchmark via `R CMD INSTALL --no-multiarch --preclean .` then a fresh
  session; `devtools::load_all()` compiles Rcpp at `-O0` and misreports C++
  speed.
