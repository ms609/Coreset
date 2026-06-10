# Red-team focus areas — MaxMin

Rotation table for the `/red-team` skill. One area per invocation. Each area pairs an algorithm's
R wrapper with its C++ core so the finder sees the `.Call`/Rcpp interface contract in one pass.

`start_tier` defaults to `sonnet` for every area — **maturity is measured, not assumed** (see SKILL.md
§Model tiers). Raise an area's `start_tier` only once cheap sweeps have stopped yielding and a higher
tier has confirmed only subtle bugs remain. Do not pre-declare maturity.

| # | Area | Files | start_tier | Key questions |
|---|------|-------|-----------|---------------|
| 1 | DropAdd (single-fidelity) | `R/dropadd.R`, `src/dropadd.cpp` | sonnet | Index bookkeeping across the R↔C++ boundary (0- vs 1-based)? Drop-then-add invariants when the candidate set empties or collapses to one point? Distance-matrix updates after a swap — stale entries? Behaviour when `n == k` or `k == 0`? |
| 2 | DropAdd (multi-fidelity) | `R/dropadd.R`, `src/dropadd_mf.cpp` | sonnet | Multi-fidelity bookkeeping vs the single-fidelity path — divergent assumptions? Fidelity-level weighting / cost accounting correct? Shared R wrapper feeding both cores — argument mismatch? (See memory: test-dropadd-mf FP flakiness — is it environment noise or a real ordering bug?) |
| 3 | GRASP | `R/grasp.R`, `src/grasp.cpp` | sonnet | Greedy-randomised construction: RCL (restricted candidate list) bounds and ties? Local-search termination — can it loop or stop early? RNG seeding/restart reproducibility? (See memory: test-grasp FP flakiness.) Empty/degenerate candidate handling? |
| 4 | Gonzalez (farthest-point) | `R/gonzalez.R` | sonnet | Farthest-first traversal correctness; first-seed choice and ties? Distance updates O(n) per step — stale minima? Degenerate inputs (duplicate points, single point, `k >= n`)? Pure-R numeric edge cases (Inf/NaN distances)? |
| 5 | Exact + maximin core | `R/exact.R`, `src/maximin.cpp`, `src/maximin_points.cpp` | sonnet | Exact/branch-and-bound search: pruning-bound correctness (can it prune the optimum)? Distance computation in the core — squared vs true distance consistency? Integer overflow on combinatorial counts? Tie-breaking vs the heuristics? |
| 6 | Seed + score | `R/seed.R`, `R/score.R` | sonnet | Seeding strategies feeding the algorithms — do they honour constraints? Score computation: is the reported max-min objective consistent with what each algorithm optimises? Off-by-one in score indexing? Empty/NA handling? |
| 7 | Test-suite health | `tests/testthat/*.R` | sonnet | Tautological or under-constrained assertions? Tests that pass on stale snapshots? FP bit-identity comparisons that are environment-fragile (see memory)? Coverage gaps for the edge cases in areas 1–6? Are claimed invariants actually checked? |

**Cross-cutting note.** The Rcpp boundary (`R/RcppExports.R`, `src/RcppExports.cpp`) is generated — don't
review it as its own area, but scrutinise the hand-written `.Call` contracts within each algorithm's area.
`R/utils.R`, `R/zzz.R`, `R/MaxMin-package.R` are small/boilerplate — fold any incidental findings into the
nearest area rather than rotating on them.
