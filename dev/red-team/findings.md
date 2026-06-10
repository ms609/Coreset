# Red-team findings — MaxMin

Non-trivial bugs and performance issues, filed **only after verification** (trivial issues are noted in
`log.md`, not here). A finding reaches a table below only if a severity-matched verifier marked it REAL
(SKILL.md §Severity-matched verification). Severity/verifier in the last column.

---

## ⭐ Headline bugs (fix-worthy)

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| FF-001 | **med-high** | FIXED 9001 | Bug | **NA/`-Inf` in distance matrix or oracle column → silent duplicate selection** | [R/farfirst.R:8-14](../../R/farfirst.R) — `.AsDistMatrix()` checks only `is.matrix`/`is.numeric`/`nrow==ncol`, **no `anyNA` guard** (unlike `.AsPointsMatrix()` at :48). `pmin.int`/`which.max` propagate NA; `which.max` returns an already-`-Inf`-masked index → a **repeated index** in the result. Reproduced: `FarFirst(matrix(c(0,5,NA,5,0,NA,NA,NA,0),3,3),3,seed=1L)` → `c(1,2,1)`. `DropAdd` matrix path silently accepts too. The `points=` path correctly errors. Realistic trigger: a tree/graph-distance oracle returning NA for disconnected pairs. **Fix: add `if (anyNA(d)) stop(...)` to `.AsDistMatrix`** (and the oracle `.DistColumn`). Surfaces from 3 angles → also T7-08 (test gap), MF-01 (NaN direct-call). | REAL (opus) |
| GRASP-01 | **HIGH** | FIXED 9001 | Bug | **Segfault at documented `alpha=1` (pure greedy) via FP-rounding empty RCL** | [src/grasp.cpp:124-132](../../src/grasp.cpp) / [R/grasp.R:72-77](../../R/grasp.R) — `thresh = gmin + alpha*(gmax-gmin)` at `alpha=1` can exceed `gmax` by ~1 ULP (seed-16: diff +2.22e-16), so the `g>=thresh` filter empties the RCL; C++ then indexes `rcl[0]` on an empty vector → **uncatchable SIGSEGV** (kills the R session). R path throws a catchable error. Live SIGSEGV confirmed: `Grasp(d,m=5,alpha=1,seed=16)` exit 139, deterministic. **Fix: empty-RCL fallback to the unique greedy-best in BOTH `grasp_construct` and `.GraspConstruct`** (preserve bit-parity). | REAL (opus), HIGH |
| GRASP-02 | **HIGH** | FIXED 9001 | Bug | No `alpha` validation in `Grasp()` → segfault for `alpha>1`, silent for `alpha<0` | [R/grasp.R:309-333](../../R/grasp.R) — every other arg is guarded; `alpha` is passed straight through. `alpha>1` empties the RCL deterministically → same `rcl[0]` SIGSEGV (confirmed `alpha=2,seed=1` exit 139); `alpha<0` silently selects from all candidates. **The empty-RCL guard (GRASP-01) plus a `stopifnot(alpha>=0, alpha<=1)` closes both.** | REAL (opus), HIGH |

---

## Round 2 — area #2 DropAdd (multi-fidelity)

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| MF-01 | low | FIXED 9001 | Bug | `DropAdd_points_cpp` accepts NaN via direct `:::` call → possible OOB write `in_S[-1]` | [src/dropadd_mf.cpp](../../src/dropadd_mf.cpp) — public `DropAdd()` is guarded by `.AsPointsMatrix` (`anyNA`), but the `[[Rcpp::export]]`ed core has no internal NaN guard; argmax degrades, `x_new` can stay `-1`. Only reachable via direct internal call. | REAL (haiku) |
| MF-03 | low | OPEN | Perf | `timeBudgetS` overshoot up to ~256 iterations at large n | [src/dropadd_mf.cpp:201](../../src/dropadd_mf.cpp) — timer checked only every `check_every=256`; at manuscript-scale n each iter is O(n·dim), so a 1 s budget can overshoot by seconds. Magnitude undocumented. | REAL (haiku) |
| MF-04 | low | OPEN | Bug | Signed-int overflow UB on `iters_done` at `max_iter=.Machine$integer.max` | [src/dropadd_mf.cpp:183](../../src/dropadd_mf.cpp) (also dropadd.cpp) — `int` counter wraps → UB; practically unreachable (billions of iters). Use `int64_t`/`size_t`. | REAL (haiku) |

## Round 3 — area #3 GRASP (see also headline GRASP-01/02)

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| GRASP-04 | low | FIXED 9001 | Bug | `.GraspConstruct` returns length-2 for `m=1` in R (`2L:m` counts *down*) | [R/grasp.R:65](../../R/grasp.R) — `for (h in 2L:m)` with m=1 is `c(2,1)`. C++ safe. Blocked for users by the public `m>=2` guard; unsafe if the helper is called directly. Fix: `seq.int(2L,m)` or an `m<2` guard. | REAL (self) |
| GRASP-05 | info | FIXED 9001 | Gap | R/C++ parity test fixes `alpha=0.8`; never exercises `alpha=1` | [test-grasp.R:26-49](../../tests/testthat/test-grasp.R) — adding `alpha=1` to the grid would have caught GRASP-01 as a test failure rather than a session crash. | REAL |

## Round 4 — area #4 FarFirst (see also headline FF-001)

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| FF-002 | low | FIXED 9001 | Bug | `seed=NA/NaN` reaches C++ as `INT_MIN` → opaque error | [R/farfirst.R:247-248](../../R/farfirst.R) — `is.integer(NA_integer_)` is TRUE so `as.integer(NA)`→NA→INT_MIN. Add `anyNA`/`is.finite` check. | REAL (haiku) |
| FF-003 | low | FIXED 9001 | Bug | `seed=integer(0)`/`c(1L,2L)` → opaque Rcpp "Expecting a single value" | [R/farfirst.R:247-248](../../R/farfirst.R) — missing `length(seed)==1` check. | REAL (haiku) |
| FF-004 | low | FIXED 9001 | Bug | `n=Inf` → spurious base-R coercion warning before package error | [R/farfirst.R:286](../../R/farfirst.R) — pre-check `is.finite(n)` before `as.integer`. | REAL (haiku) |
| FF-005 | info | OPEN | Design | Ensemble returns bare `seq_len(nPts)` (no `strategy_results`/`winning_strategy`) when `n>=nPts` | [R/seed.R:279,378](../../R/seed.R) — inconsistent with documented "with attributes" return. (= F-603.) | REAL (haiku) |
| FF-006 | info | OPEN | Design | Asymmetric distance matrix silently accepted | [R/farfirst.R:8-14](../../R/farfirst.R) — `nrow==ncol` only, no symmetry check. O(N²) to check; if intentional, document. | REAL (haiku) |

## Round 5 — area #5 Exact + maximin core  (seam ran DRY at sonnet → next visit escalates to opus)

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| E5-01 | low | FIXED 9001 | Bug | `MaximinFrom_cpp`/`MaximinFromPoints_cpp` lack an internal `1<=n<=nPts` check | [src/maximin.cpp:13](../../src/maximin.cpp), [src/maximin_points.cpp:92](../../src/maximin_points.cpp) — reachable only via direct `:::`; public wrappers guard. Defensive `Rcpp::stop`. | REAL (haiku) |

## Round 6 — area #6 Seed + score

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| F-601 | low-med | FIXED 9001 | Bug | `match.arg(several.ok=TRUE)` silently drops misspelled anchor names | [R/seed.R:273-277,372-376](../../R/seed.R) — reachable via `FarFirst(seed=c("peripheral","anti-medoid"))` (FarFirst doesn't validate a multi-element seed). User silently gets a smaller ensemble; `winning_strategy` gives no hint. Reproduced. Add an explicit `setdiff` check. | REAL (opus) |
| F-602 | low | FIXED 9001 | Doc | `winning_strategy` reports only the first strategy when all `t_k` are NA (n==1) | [R/seed.R:97-99](../../R/seed.R) — doc claims "all tied-best strategies". Low stakes (n==1 is trivial). | REAL (haiku) |
| F-604 | low | FIXED 9001 | Bug | `MinDist`/`.SubsetScore` with duplicate `idx` → returns 0 (self-distance survives `diag<-Inf`) | [R/score.R:29](../../R/score.R), [R/farfirst.R:118](../../R/farfirst.R) — no guard; a repeated index makes any selection score worst-possible. Caller-contract issue. | REAL (haiku) |
| F-605 | low | FIXED 9001 | Bug | NA in `idx`: matrix path returns NA silently; coordinate path errors — inconsistent | [R/farfirst.R:118](../../R/farfirst.R) — add `if (anyNA(idx)) stop(...)` to both. | REAL (haiku) |

## Round 7 — area #7 Test-suite health

| ID | Severity | Status | Type | Title | Location & detail | Verifier |
|----|----------|--------|------|-------|-------------------|----------|
| T7-08 | high (gap) | FIXED 9001 | Gap | No test for NA in a distance matrix (= FF-001 root) | [tests/testthat/](../../tests/testthat) + `.AsDistMatrix` — add `expect_error(DropAdd(na_mat,m=3),"NA")` / `FarFirst(na_mat,3)` once the guard lands. | REAL (haiku) |
| T7-01 | med | OPEN | Weak | Path-relink `>= max(endpoints)` assertion holds by construction | [test-grasp.R:106](../../tests/testthat/test-grasp.R) — `bestZ` initialised to `max(zX,zY)` and only increases. Use `>` on a geometry where an intermediate improves. | REAL (haiku) |
| T7-04 | med | OPEN | Weak | DropAdd seed-reproducibility test is vacuous (no RNG draws) | [test-dropadd.R:190-193](../../tests/testthat/test-dropadd.R), test-dropadd-mf.R:149-151 — DropAdd is deterministic; same-seed identity is trivial. (The dead `seed` arg is being removed by the API task.) | REAL (haiku) |
| T7-06 | med | OPEN | Weak | `expect_identical` on float `score`/`secondary` (FP-fragile; root of the known flakiness) | [test-dropadd.R:339-355](../../tests/testthat/test-dropadd.R) — split: `expect_identical` on integer indices, `expect_equal(tolerance=)` on floats. | REAL (haiku) |
| T7-07 | med | OPEN | Weak | `expect_identical` on float `score` (FP-fragile) | [test-grasp.R:45](../../tests/testthat/test-grasp.R) — same fix as T7-06. | REAL (haiku) |
| T7-09 | med | OPEN | Gap | `secondary` attribute never checked against the brute-force upper-triangle sum | [test-dropadd.R](../../tests/testthat/test-dropadd.R) — parity tests catch cross-path divergence, not a shared formula bug. | REAL (haiku) |
| T7-10 | med | OPEN | Gap | ExactMaxMin `proven=FALSE` (budget-expiry) branch untested | [test-exact.R](../../tests/testthat/test-exact.R) — add a `timeBudgetS=1e-6` test asserting `!proven`. | REAL (haiku) |
| T7-12 | med | OPEN | Weak | farfirst test strips attributes before compare, hiding `strategy_results` divergence | [test-farfirst.R:34-35](../../tests/testthat/test-farfirst.R) | REAL (haiku) |
| T7-03 | low | OPEN | Weak | Exact `m=n` → `1:8` is structurally forced | [test-exact.R:91](../../tests/testthat/test-exact.R) — assert `objective` instead. | REAL (haiku) |
| T7-14 | low | OPEN | Gap | FarFirst `n>N` tail order (positions 6..N) untested | [test-farfirst.R](../../tests/testthat/test-farfirst.R) | REAL (haiku) |
| T7-15 | low | OPEN | Gap | Grasp `eliteSize=1` (PR skipped, `pr_calls=0`) untested | [test-grasp.R](../../tests/testthat/test-grasp.R) | REAL (haiku) |
| T7-16 | low | OPEN | Gap | `.GraspPathRelink` identical-input early return untested | [R/grasp.R:155-157](../../R/grasp.R) | REAL (haiku) |

---

### Trivial / cosmetic (noted, not fix-worthy on their own — fold into the relevant fix)
~~GRASP-03 dead var `candMask` (grasp.R:66)~~ (FIXED 9001) · ~~F-606 cosmetic guard `n<1L` vs `n==0L` (seed.R:280/379)~~ (FIXED 9001) · FF-007 oracle default-seed doc · FF-008 non-integer `n` truncation · MF-02 want_trace coverage artifact (← DA-01) · E5-02 Recover() first-m (valid lower bound, non-bug) · T7-02/T7-05 redundant assertions · T7-11/T7-13/T7-17 minor test nits.
