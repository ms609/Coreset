# Agent notes for `MaxMin`

`MaxMin` is a small, dependency-light (`Imports: Rcpp, stats`) solver library
for the Max-Min Diversity Problem (MMDP / discrete p-dispersion). It was
extracted from the study package `FurthestPoint` so that CRAN packages (notably
`TreeSearch`, via `TreeSearch::WideSample()`) can depend on the solvers without
pulling in manuscript machinery. **Keep it CRAN-clean and dependency-light** —
do not add `Imports` beyond `Rcpp`/`stats`; `highs` (exact solver) stays in
`Suggests` behind a `requireNamespace()` guard.

## Scope

In scope: the MMDP algorithms themselves — Gonzalez farthest-first (matrix /
coordinate / column-oracle paths) and its seeding strategies, DropAdd tabu
search, the exact node-packing solver, the polish local search, and the
`TkScore` objective. Out of scope: benchmark harnesses, competitor comparisons,
synthetic/real test-case constructors, manuscript figures — those live in
`FurthestPoint`, which `Imports: MaxMin`.

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

## The correctness contract

The matrix, coordinate, and column-oracle paths must return **bit-identical**
selections on the same data (the kernels reproduce `stats::dist()`'s Euclidean
bits exactly). `FurthestPoint` reproduces published manuscript numbers against
this package, so changing any solver's tie-breaking or seeding silently is a
breaking change — guard it with the reproduction tests there.

## Not `maximin`

Unrelated to the CRAN package `maximin` (space-filling designs). Keep the "See
also" cross-reference in the package documentation.
