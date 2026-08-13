# Profiling log — MaxMin

## Round 1 — 2026-06-10 — Area 1: FarFirst greedy + ensemble

**Question (user):** Is `MaximinFrom_cpp` the right C++ entry point, or should the
R "wrangling" move into C++ so the R↔C++ barrier is crossed once (with speedup)?

**Driver:** `drivers/farfirst-default.R` — default coordinate path
(`seed="random_furthest"`, 3 starts), phase breakdown (seed-find / greedy pass /
MinDist scoring) vs n/N, N=6000, dim ∈ {2,10}, n/N ∈ {0.01,0.05,0.20,0.50}.
**Verification:** `drivers/verify.R` (1620-case correctness + timing battery,
run against the pre-edit baseline lib and the patched lib) + `drivers/compare.R`.

**profvis/breakdown verdict:**
- Seed-find + R wrangling: **< 0.2 ms total, independent of n** → barrier-crossing
  is NOT the cost. The "fold into C++" premise is refuted (finding T-003).
- `MinDist` re-scoring scales as ~n/N of the pass: 3% (n/N=0.01) → 78% (n/N=0.5).
  Redundant — the greedy pass already computes T_k (finding T-001).
- The pass itself is the irreducible core; its inner loop is sqrt-bound at low
  dim (finding T-002).

**Implemented (T-001, T-002):**
- `src/maximin.cpp`, `src/maximin_points.cpp`: track running min insertion
  distance, attach as `t_k` attribute (free); points kernel works in squared
  space, single final sqrt. New `EuclidColSq` helper.
- `R/gonzalez.R`: `.StripScore()`; bare single-pass strips `t_k`. `R/seed.R`:
  both ensemble drivers read `attr(idx,"t_k")` instead of `MinDist`.

**Result (verified, end-to-end FarFirst, 3 starts, N=6000):**
```
 dim    n  n/N | points old->new        | matrix old->new
   2   60 0.01 |  4.1 -> 3.5  (1.17x)    |  3.0 -> 2.7  (1.14x)
   2  300 0.05 | 17.1 ->14.0  (1.23x)    | 14.4 ->10.1  (1.43x)
   2 1200 0.20 | 81.2 ->52.7  (1.54x)    | 88.8 ->36.9  (2.40x)
   2 3000 0.50 |280.4 ->128.5 (2.18x)    |362.9 ->90.0  (4.03x)
  10   60 0.01 |  9.5 -> 7.8  (1.22x)    |  3.3 -> 2.7  (1.24x)
  10 3000 0.50 |784.4 ->322.1 (2.44x)    |368.7 ->92.0  (4.01x)
```
Correctness: cross-path identical 1620/1620 (NEW & OLD); OLD==NEW selection
1620/1620 both paths; free-T_k vs MinDist max abs err 0.000e+00. Full suite:
the only failures (test-dropadd ×4, test-gpr ×17) are pre-existing FP/RNG
bit-identity tests in solvers untouched by this change (confirmed identical on
the baseline lib).

**Open:** T-004 (BLAS/Gram gemv for the column) — bigger lever on the pass,
trades the cross-path guarantee; deferred to a future round.

Status: Area 1 → OPTIMISED.

---

## Round 2 — 2026-06-10 — Area 2: DropAdd tabu search

**Question:** Where does DropAdd's time go — across n, m (small vs large), and
matrix vs points — and what is the per-iteration cost structure?

**Driver:** `drivers/dropadd.R` — matrix path n∈{2000,4000} m∈{10,n/2};
points path n∈{20000} dim∈{2,10} m∈{10,n/2}, plus n=4000 dim=2 for cross-path
comparison. Fixed `maxIter=1500, plateau=.Machine$integer.max` (exactly 1500 iters
per case). Bare worst case: n=20000 dim=10 m=10000 ≈ 2.6 s/rep.

**Timing table (median of 5–10 reps, 1500 iters):**

| path   |     n | dim |     m | ms_total | ms/iter |
|--------|-------|-----|-------|----------|---------|
| matrix |  2000 |  —  |    10 |     35.0 |  0.0233 |
| matrix |  2000 |  —  |  1000 |     70.0 |  0.0467 |
| matrix |  4000 |  —  |    10 |    160.0 |  0.1067 |
| matrix |  4000 |  —  |  2000 |    250.0 |  0.1667 |
| points |  4000 |   2 |    10 |     60.0 |  0.0400 |
| points |  4000 |   2 |  2000 |    110.0 |  0.0733 |
| points | 20000 |   2 |    10 |    350.0 |  0.2333 |
| points | 20000 |   2 | 10000 |   1460.0 |  0.9733 |
| points | 20000 |  10 |    10 |    830.0 |  0.5533 |
| points | 20000 |  10 | 10000 |   2570.0 |  1.7133 |

No case exceeded 8 s (worst single rep ≈ 2.6 s).

**CORRECTED analysis: construction vs search subtraction (`maxIter=0`):**

Construction costs confirmed by timing with `maxIter=0L`:

| path   |     n | dim |     m | construct_ms | search_ms | ms/iter (search) |
|--------|-------|-----|-------|-------------|-----------|-----------------|
| matrix |  2000 |  —  |    10 |        20.0 |      30.0 |          0.0200 |
| matrix |  2000 |  —  |  1000 |        30.0 |      40.0 |          0.0267 |
| matrix |  4000 |  —  |    10 |       115.0 |      45.0 |          0.0300 |
| matrix |  4000 |  —  |  2000 |       160.0 |      90.0 |          0.0600 |
| points |  4000 |   2 |    10 |         0.0 |      65.0 |          0.0433 |
| points |  4000 |   2 |  2000 |        30.0 |      80.0 |          0.0533 |
| points | 20000 |   2 |    10 |         0.0 |     350.0 |          0.2333 |
| points | 20000 |   2 | 10000 |      1060.0 |     390.0 |          0.2600 |
| points | 20000 |  10 |    10 |         0.0 |     830.0 |          0.5533 |
| points | 20000 |  10 | 10000 |      1880.0 |     680.0 |          0.4533 |

**Scaling (search only, per-iter costs):**
- n scaling, matrix m=10: n=2000→4000 per-iter 0.0200→0.0300ms = ×1.5 — consistent
  with O(n) search (two full column scans per iter). NOT ×4.6 — the total-time
  super-linearity was the O(n²) construction seed dominating the ms_total.
- m effect, matrix search: n=4000 m=10→2000: 0.0300→0.0600ms = ×2.0 — reflects the
  O(m) recompute scan over S and the O(m) objective re-evaluation each iter.
- n scaling, points dim=2 m=10: n=4000→20000 (5x): 0.0433→0.2333ms = ×5.4 — clean
  O(n) (two EuclidCol column passes per iter).
- m effect, points search: n=20000 dim=2 m=10→10000: 0.2333→0.2600ms = ×1.1 (flat!).
  n=20000 dim=10 m=10→10000: 0.5533→0.4533ms = flat. **The recompute branch is NOT
  costly at these sizes** — the search per-iter cost is dominated by the two O(n·dim)
  column passes regardless of m.
- dim multiplier, points n=20000 small m: dim=2→10 = 0.2333→0.5533ms = ×2.4
  (not ×5; SIMD vectorises the inner dim loop, especially at dim=2 which gets one
  SSE or scalar pair, while dim=10 gets at least 2 SIMD passes).
- **Construction dominates at large m + large n (points):** n=20000 m=10000
  construction is 1060–1880 ms vs 390–680 ms for 1500 search iters. The O(m·n·dim)
  greedy construction greedy loop is the bottleneck for large-m matrix-free calls.
- **Matrix construction is O(n²):** the seed computes n full row sums (`dropadd.cpp`
  lines 50–58); at n=2000 this is 20 ms, n=4000 is 115 ms — ratio 5.5x vs expected
  4x (memory bandwidth compounds the n² work at large matrices).

**profvis verdict (matrix n=4000 m=2000, 5 reps; points n=20000 dim=2 m=10000, 3 reps):**
- Matrix: DropAdd_cpp (C++) = 98.8% of all profvis intervals. No R label >2%.
  (The DropAdd R frame appears in same intervals as its .Call leaf — this is normal.)
- Points dim=2: DropAdd_points_cpp (C++) = 100%. No R overhead.
- Points dim=10: DropAdd_points_cpp (C++) = 99.8%. No R overhead.
- Conclusion: ALL cost is inside the C++ kernels; no R hotspot is actionable.

**Finding (T-005) — REVISED:**
The actionable targets are both in the construction phase of `DropAdd_cpp` /
`DropAdd_points_cpp`, not the tabu search loop:
1. Matrix path `dropadd.cpp` lines 50–58: O(n²) seed computation (double loop over
   all (i,j) pairs for row-sums). At n=4000, construction is 115 ms vs 90 ms for
   1500 search iters — construction ≈ 56% of total call cost. The O(n²) seed can be
   replaced by the same O(n·dim) anti_centroid-distance proxy already used in
   `dropadd_mf.cpp`.
2. Points path `dropadd_mf.cpp` lines 119–165: O(m·n·dim) greedy construction loop;
   at n=20000 m=10000, construction is 1060–1880 ms vs 390–680 ms for 1500 iters.
   A fast O(n·dim) construction would cut total call time in half at large m.
VTune is still needed for line-level attribution inside the search loop, but
construction is the more impactful target confirmed at R level.

Status: Area 2 → PROFILED.

last_focus: 3

---

## Round 3 — 2026-06-10 — Area 3: Grasp refinement

**Question:** Where does Grasp's time go — construction (phase A), GRASP
iterations + local search (phase B), or path relinking (phase C) — at small m
and large m?

**Driver:** `drivers/grasp.R` — dense-matrix only; n∈{200,500}; m∈{10,50,100}
at n=200; m=10 only at n=500 (large-m n=500 >> 8 s, shrunk/omitted). eliteSize
varied {3,5,8} at n=200/m=100 to isolate PR contribution (PR ∝ eliteSize²).
plateau=15, seed=1L throughout. Timing: median of 5 calls per case. n=500
large-m cases shrunk per ≤8 s constraint (noted in driver header).

**Timing table (median of 5 calls):**

| label             |   n |   m | eliteSize | plateau | iters | pr_calls | ms_median |
|-------------------|-----|-----|-----------|---------|-------|----------|-----------|
| n200_m10_es5      | 200 |  10 |         5 |      15 |    28 |       20 |       5.5 |
| n200_m50_es5      | 200 |  50 |         5 |      15 |    15 |       20 |     440.0 |
| n200_m100_es3     | 200 | 100 |         3 |      15 |    18 |        6 |    1010.0 |
| n200_m100_es5     | 200 | 100 |         5 |      15 |    16 |       20 |    2200.0 |
| n200_m100_es8     | 200 | 100 |         8 |      15 |    15 |       56 |    5450.0 |
| n500_m10_es5      | 500 |  10 |         5 |      15 |    30 |       20 |      14.5 |

(n=200 m=10 timed via block of 20 calls; OS timer resolution ~16 ms masks individual calls.)

**Scaling: m (small vs large) at fixed n=200, eliteSize=5:**
- m=50→100 (×2 in m): time 440→2200 ms = ×5.0, empirical m-exponent ≈ 2.3
- Consistent with O(m²) local-search scoring plus O(m²) PR candidate scoring,
  both growing faster than linear but slower than cubic.

**PR isolation (eliteSize variation at n=200, m=100, plateau=15):**
- eliteSize 3→5: +14 PR calls → +1190 ms = **85.0 ms/PR call**
- eliteSize 5→8: +36 PR calls → +3250 ms = **90.3 ms/PR call**
- Mean PR call cost: **~88–96 ms** at m=100 (consistent ≈ linear in pr_calls)
- At m=10: all 20 PR calls total <1 ms (PR negligible at small m)

**Direct phase discriminator (maxIter=0 at n=200, m=100, eliteSize=5):**

| condition     | ms_median | iters | pr_calls | meaning            |
|---------------|-----------|-------|----------|--------------------|
| maxIter=0     |      2070 |     0 |       20 | phase A + phase C  |
| full run      |      2160 |    16 |       20 | phase A+B+C        |
| eliteSize=1   |        10 |     0 |        0 | phase A only       |

Phase decomposition (eliteSize=5, m=100):
- Phase A (5 elite builds, ~10 ms each): **50 ms (2.3%)**
- Phase B (GRASP loop, ~16 iters × ~5.6 ms): **90 ms (4.1%)**
- Phase C (20 PR calls, ~100 ms/call): **2020 ms (93.6%)**

Note: each PR call in phase C includes the PR walk + 1 local-search pass on
the result. Both the walk (m×m sub-matrix min per step) and the LS (m×m sub-
matrix min per swap candidate) use the same O(m²) re-score primitive.

**profvis verdict (n=200, m=100, eliteSize=5, plateau=15, 3 reps ≈ 6.6 s):**
- 418/418 time-point samples have `Grasp_cpp` at depth 2 (inside the .Call) →
  **100% of wall time inside the C++ kernel**
- R wrapper (set.seed + .AsDistMatrix + arg checks): 0 samples → unmeasurably fast
- No R label > 2%. All actionable cost is in `Grasp_cpp`.

**Dominant phase by regime:**
- **Small m (m=10):** All phases cheap (PR negligible). Total ~5.5 ms/call; PR
  share <1 ms, construction + LS ~10× faster (m=10 → m²=100 vs m=100 → m²=10000).
- **Large m (m=100):** Phase C (PR) = **93.6%** of total. Phase B = 4.1%. Phase A = 2.3%.
  Each PR call re-scores m×m sub-matrices at every step of the PR walk and after
  the LS on the PR result — the shared O(m²) re-score primitive is the hot kernel.

**Finding (T-006) — REVISED:**
`[Optimise]` The shared sub-matrix min-rescore primitive in `Grasp_cpp` drives
93.6% of wall time at large m (phase C) and 4.1% (phase B), both via the same
O(m²) `min(d[cand,cand])` re-score. Replacing this with an incremental swap update
(drop one row+col, add one row+col, maintain running row-min and global min → O(m)
per candidate) would reduce both the PR walk and the local-search swap scan. The
fix has two call sites (PR step loop and LS swap scan); both benefit. VTune line-
level on `Grasp_cpp` to confirm the sub-matrix min scan is the dominant hot line.

Status: Area 3 → PROFILED (then OPTIMISED in Round 4).

---

## Round 4 — 2026-06-10 — Implementation of T-004, T-005, T-006

Implemented the three findings from Rounds 1–3, one per file (no collision):
`src/dropadd.cpp` (T-005a), `src/grasp.cpp` (T-006), `src/maximin_points.cpp`
(T-004). All three are **bit-identical** to the prior build by construction, so
verification = "results match the baseline lib exactly". Drivers
`drivers/verify-r2.R` (per-lib battery) + `drivers/compare-r2.R`.

**Correctness — 416/416 bit-identical vs the prior (T-001/T-002) build:**
FarFirst 128/128, DropAdd matrix 48/48, DropAdd points 48/48, Grasp 192/192.

**Timing — old → new:**
```
T-005 DropAdd matrix construction (maxIter=0):
  n=4000   116.1 ->  10.3 ms   (11.3x)
  n=6000   293.4 ->  22.4 ms   (13.1x)
T-006 Grasp (n=200):
  m=50     315.5 ->  22.2 ms   (14.2x)
  m=100   1328.8 ->  85.0 ms   (15.6x)
T-004 FarFirst points single pass (column reorder):
  d2  n300/1500/3000   1.06–1.08x      (consistent small low-dim win)
  d10 n300/1500/3000   0.99–1.01x      (noise — pass at limit)
```

**Implementations:**
- T-005a `dropadd.cpp`: matrix seed row-sum accumulated column-sequentially
  (`for j: for i: rs[i]+=col[i]`) — same summation order, sequential access.
- T-006 `grasp.cpp`: PR-walk and local-search candidate scoring hoist the
  per-drop `min-pair-within(rem)` out of the add-loop; each add is O(m) via the
  new `min_to_set`. LS pair-count (`count_pairs_le`) deferred to nd ≥ best only.
  `eval_cand` removed (dead). `min` is exact ⇒ bit-identical.
- T-004 `maximin_points.cpp`: new `EuclidColSqInto` fills the squared-distance
  column dimension-by-dimension over contiguous columns into a reused buffer;
  `EuclidColSq` (per-point) removed. Bit-identical; marginal — area AT-LIMIT.

**Full suite vs final install:** test-gonzalez 88/0, test-seed-score 45/0,
test-exact 51/0, test-dropadd-mf 39/0 all green. The only failures are the
pre-existing FP "C++ vs R-reference bit-identical" tests: test-dropadd
**fluctuates 4↔5 run-to-run** (FP-environment-sensitive across random instances,
confirmed by repeated runs of the SAME build) and test-gpr 17 — neither caused
by these changes (results are 416/416 bit-identical to baseline).

Status: Area 2 → OPTIMISED, Area 3 → OPTIMISED, Area 1 pass → AT-LIMIT.

last_focus: 3

Status: Area 3 → PROFILED.

## Round (T-007) — 2026-06-10 — Area 3: Grasp base_z hoist (from salvaged audit)

**Origin:** an independent peak-optimality audit (subagent, killed mid-run by
token depletion before its verdict) surfaced one un-chased lead: in
`grasp_path_relink`, `base_z = objective_of(pk \ {di_})` is recomputed O(m²)
**per drop candidate** — T-006 had only hoisted it out of the inner *add* loop.
Confirmed by reading the kernel; the same pattern lives in `grasp_local_search`.

**Fix:** new `min_edge_witness(sel, &wa, &wb)` returns the global min pairwise
edge of `sel` and the two vertices realising it (O(m²), once per step). For any
drop `di_ ∉ {wa,wb}` the witness edge survives in `pk\{di_}`, so the post-drop
min is unchanged → `base_z = gmin` in O(1); only the ≤2 witness vertices rescore.
Applied at `grasp.cpp:245` (PR walk) and `:159` (LS, reusing the di-loop's min).
Per-step base_z: O(|drop_cands|·m²) → O(m²). Size guard kept first.

**Why bit-identical:** `min` re-selects the identical surviving `D()` value (same
matrix cell; `sel`/`rem` ascending, witness read at `a<b`). No arithmetic.

**Verified (verify-r2.R old vs new):**
```
correctness  416/416 identical   (gp 192/192; ff/da untouched)
grasp_n200_m50    21.92 -> 16.35 ms   (1.34x)
grasp_n200_m100   81.90 -> 39.52 ms   (2.07x)
all FarFirst/DropAdd timing flat (<=1.03x)
```
Win grows with m (PR phase = 93.6% of wall); analytic worst case r=m ≈ 2.3–2.5×
on the PR walk. Comfortably above the 5% floor.

Status: Area 3 → OPTIMISED (T-007). Rotation current: all areas OPTIMISED/AT-LIMIT.

last_focus: 3

---

## Round 5 — 2026-06-11 — Area 4: ExactMaxMin node-packing (user-requested)

**Question (user):** the manuscript's exact solver maxes out at n=240 (Fig. 1).
Profile + optimise `ExactMaxMin` until no >10% gain remains, so the six 342–990
datasets can be solved exactly within a 3-day HPC budget.

**Driver:** `drivers/exact.R` (penguins n=342 k=10, the smallest unsolved target);
`drivers/exact-prof.R` (Rprof line+memory); `drivers/exact-probe.R` (per-probe cost
instrumentation); `drivers/exact-characterize.R` (heuristic-vs-exact hit rate over
14 ground-truth cases × k∈{2,4,6,10}); `drivers/capture_old.R` + `verify_new.R`
(OLD-vs-NEW regression + RNG contract). ExactMaxMin is pure R over the `highs`
backend, so this was a profvis/Rprof (R-triage) round, no VTune.

**Profile verdict (penguins n=342 k=10):**
- `solver_run` (highs MILP, C) = **86%** self time — the IP solves dominate.
- `which(d < lambda & upper.tri(d), arr.ind=TRUE)` = 7.4% self, **933 MB** —
  two fresh n×n logical matrices *per probe*.
- dense `matrix(0, nEdge, n)` for the packing matrix = 4.4% self, 176 MB at
  n=342 → **~3.9 GB at n=990**: the actual scaling wall (the n≈240 cap).
- Per-probe instrumentation: ~16 evenly-costed probes (~0.4 s each), ~half of
  them splitting near-identical thresholds (2.2125 / 2.2140 / 2.2143 …) around
  the optimum — each a full IP solve.

**Levers tried (measured before filing — see findings T-008):**
1. **Feasibility reformulation** (constant obj + `sum(x)≥m` cut, stop at first
   feasible). **REJECTED**: 0.73–0.89× (slower). The cardinality cut weakens the
   LP relaxation and removes the objective guidance highs prunes with;
   `objective_bound` cutoff had no effect under `maximum=TRUE`.
2. **Sparse packing matrix** (`Matrix::sparseMatrix`, 2 nnz/edge). ~1.15× on its
   own AND removes the GB-per-probe memory wall → the enabler for n≥683.
3. **Heuristic warm-start gallop.** Seed a provably-achievable lower bound (best
   of 8 `Grasp` restarts + `DropAdd`), gallop up to the first infeasible
   threshold, bisect the small bracket. Grasp attains the optimum where DropAdd
   does not; on the ground-truth grid the heuristic == exact in **55/56**
   (case,k) → a single infeasibility solve certifies it (k=2: zero solves).

**Correctness is independent of seed quality** — the warm start only sets the
starting lower bound; the search proves the true optimum regardless (the one
heuristic miss, tc6 k=10, still returned the exact optimum, in 9 probes vs 1).
Verified bit-identical: 56/56 proven optima match the dense solver; brute-force
oracle (test-exact.R) + full suite 598 PASS / 0 FAIL.

**Result (verified, `verify_new.R`):** grid total OLD 222.6 s → NEW 11.3 s
(**19.6×**), per-case median 16.3×, max ~1300× (k=2). RNG contract: `objective`
seed-independent, selection reproducible under `set.seed()` (matches `Grasp()`);
ExactMaxMin advances the session RNG like other randomised routines (test #4
updated from the old "two calls give identical indices" to this contract, per
user steer that the exact *objective* — not the witness — is what is determined).

**Implementation:** `R/exact.R` (sparse `.MaxISVerdict`, `.ExactWarmStart`,
gallop+bisect main, new `warmStart=` arg); `Matrix` → Suggests; NEWS + roxygen.

Status: Area 4 → OPTIMISED (T-008). Headline 19.6× ≫ 10% floor; further levers
(clique cuts, warm-started simplex across probes) are speculative and unmeasured —
not pursued. Rotation current: all areas OPTIMISED/AT-LIMIT.

last_focus: 4

---

## Round 6 — 2026-06-11 — Area 5: KCentre CDSh kernel (user-requested, new code)

**Context:** new discrete k-centre solvers added to the package (`KCentre`/CDSh,
`ExactKCentre`, `KCentreRadius`). CDSh beats the Gonzalez 2-approximation that
`FarFirst` gives for this objective by 17–45% on real data (matches the proven
optimum at k=5/8 on an iris subset). User asked to `/profile` the new code; this
round profiles the CDSh kernel `KCentreCDSh_cpp`.

**Driver:** `drivers/kcentre.R` — clustered Gaussian n=2004, dim=10, 12 clusters
(representative threshold-graph density), k=20, 4 reps (~5 s bare). Built against
the optimised `-O2` install at `…/Temp/maxmin-optlib` (default-lib MaxMin.dll was
locked by an active long-running R process; do NOT install to the default lib).
**profvis:** `drivers/kcentre-profvis.R`. **Verify:** `drivers/kcentre-verify.R`
(OLD `optlib` vs NEW `optlib2`) + `drivers/kcentre-compare.R`.

**profvis verdict (n=2004 k=20):**
- `KCentreCDSh_cpp` (C++ kernel) = 233 samples (~71%) — the dominant cost.
- `.KCentreCandidates` (R) = ~91 samples (~28%): `unique.default` 40 + `sort` 29 +
  `order` 24 + `upper.tri` 17 — i.e. `sort(unique(d[upper.tri(d)]))` allocates an
  n×n logical mask and churns ~n²/2 values through R's unique/order.
- `.KCentrePeripheralSeeds` (FarFirst seed) = 3, coercion = ~15: negligible.

**Two stacked levers (T-010), both bit-identical (symmetric, column-major d):**
1. **Kernel cache reorder.** The score-init (`score[i]=#{j:d(i,j)≤r}`) and the
   domination-decrement inner loop read *rows* of a column-major matrix (stride-n,
   a cache miss per access — the T-005 pattern). Since `d(i,j)=d(j,i)`, read them
   through the column pointer `P+i*n` (stride-1). Bit-identical.
2. **Candidate enumeration → C++.** `KCentreCandidates_cpp` extracts the upper
   triangle column-sequentially then `std::sort`+`std::unique`, dropping the R
   logical mask and unique/order overhead. Same sorted distinct set.

**Result (verified, `kcentre-verify.R` OLD vs NEW, n=2004 k=20, 6 reps):**
```
per-call  1275 ms -> 307 ms   (4.16x)
centres + radius bit-identical at k = 5, 20, 50  (kcentre-compare.R)
```
Correctness on the optimised source: `test_local(filter=kcentre)` 70/70 green
(KCentreRadius vs brute force; CDSh ≤ Gonzalez; ExactKCentre == brute-force
optimum; edge cases). Full suite unaffected (kcentre is new; existing solvers
untouched).

**Open (T-011):** `ExactKCentre` not yet profiled — by analogy to T-008 (its dual
node-packing solver) the cost will be the `highs` covering-IP solves, with the
per-probe `which(d <= r, arr.ind=TRUE)` (n×n logical) the actionable R lever.
Next rotation target for area 5.

Status: Area 5 CDSh → OPTIMISED (T-010); ExactKCentre → PENDING (T-011).

last_focus: 5

## Round (targeted) — 2026-06-29 — MaxEntropy (maxdet) selector  [user: /profile MaxEntropy]

**Task:** profile `MaxEntropy()` (R/maxentropy.R + src/maxentropy.cpp) and optimise
until the best per-round gain < 20%.

**Triage (step-timing, not VTune):** the only O(n³) costs are the R-side
`eigen()` decomposition(s) and the clip-repair reconstruction; the C++ greedy is
O(n²k) and negligible at large n. VTune unnecessary — the hotspot is LAPACK
called from R, located by an internal-step breakdown
(`dev/profiling/drivers/maxentropy.R`).

**Round 1 finding (T-012, APPLIED, +36.3%):** two stacked levers in the new
`.MaxEntropyPrepare()` — (a) compute the symmetric `eigen()` ONCE instead of twice
(values-only for negMass + full for repair, on the same matrix); (b) reconstruct
the clip repair with `tcrossprod(V·√λ)` (dsyrk) instead of `V·(Λ·Vᵀ)` (dgemm),
3.6× faster, identical to 1.78e-15. n=1200 prep 3.17 s → 2.02 s.
Verified safe: test_local 180/180; greedy↔dev-prototype parity 12/12 cases; exact
10/11 identical + d3 equal-logdet tie (pre-existing ME-003 scorer tie-break, not
introduced by this change).

**Round 2 (T-013, AT-LIMIT, < 20% → STOP):** after T-012 the remaining cost is 88%
the LAPACK symmetric eigendecomposition (1.77 s of 2.02 s), required by the
nearest-PSD repair and algorithmically irreducible without approximating the repair
(which would break selection parity). Best non-eigen win < 12%. Stopping criterion
(< 20%) reached. Area marked AT-LIMIT.

**Cleanup:** no VTune result dirs created. last_focus left at 5 (this was a targeted
run; the ExactKCentre T-011 rotation slot is unchanged).

---

## Round — 2026-08-12 — Area 3: Grasp refinement (round 3)

**Question (user):** before re-running the manuscript's canonical grid, is `Grasp()`
as optimised as it can be? A targeted run, not a search for the next general
hotspot.

**Shape:** taken from the canonical pipeline — GRASP runs coreset-thinned at
`maxCandidates = 2000`, k ∈ {10, 100}, over a `plateau` ladder (canonical ≈ 512).
**Drivers:** `drivers/grasp-battery.R` (1458-cell bit-identity battery over
n/dim/k/eliteSize/alpha/plateau/seed, comparing indices, objective, `iters` and
`pr_calls`, and re-checking R↔C++ parity on the small cells) and
`drivers/grasp-timing.R` (clean-build timing, best-of-3).

**Method.** T-006/T-007 had already reduced both swap scans to O(m) inner terms, so
a sampling profiler could no longer separate them. Instead a throwaway instrumented
build counted inner operations directly — call and op counts for the min-to-set scan
and the tie-break pair count, plus critical-position counts.

**The finding that shaped the round:** the two remaining terms' shares *swap* with
`plateau`. At low `plateau` phase C (path relinking) is ~90% of wall time; at high
`plateau` — canonical included — phase B's local search is ~93%. Optimising either
alone would have looked good on one half of the ladder and done nothing for the
other. `crit`/pass is exactly 2.00, which is what kills the most obvious lever
(T-015).

**Round 3 finding (T-014, APPLIED, 2.03–3.79× at k=100):** (a) hoist the
add-to-`pk` minima in the PR walk behind a two-smallest summary, dropping a factor
of `|drop_cands|`; (b) split the LS extended-improvement pair count into its rem×rem
and s×rem halves and memoise the former per (ci, thr) — 72% hit rate. Verified
1458/1458 bit-identical including `iters` and `pr_calls`; R↔C++ parity 324/324; full
suite 4893 pass / 0 fail; all 7 new C++ branches exercised.

**Rejected en route (T-015):** the same hoist inside the local search. Because
`crit` is exactly 2 it can save ~2× in op count at best, and `near_two` costs more
per element than the branchless min it would replace — net +1.17× at k=100 but
**0.64× at k=10**, a regression on a canonical cell. A branchless row-buffer rewrite
meant to fix that was worse again (16.05 s vs 11.13 s at plateau 256): the strided
reads it targeted are already cheap, because consecutive candidates share cache
lines down each column. Both rejected variants were themselves bit-identical, so
this was a pure cost decision. Recorded so it is not re-attempted.

**Note on timings.** These are local *relative* ratios (clean build vs clean build,
best-of-3), used only to choose between implementations. Canonical wall-clock
figures for the manuscript belong on Hamilton, not this Windows box.

**Cleanup:** scratch libraries and instrumented builds live under the session
scratchpad, outside the repo; no instrumentation ships. last_focus unchanged (this
was a targeted run, not a rotation slot).

---

## Round — 2026-08-13 — Area 3: Grasp round 4 (bit-identity lifted; kept anyway)

**Question (user, via ROUND4_PLAN.md):** improve GRASP as much as possible with
the bit-identity constraint explicitly lifted; quality per unit time is the
metric; the designated lever is an incrementally-maintained nearest-selected
summary (its §4), whose cost was assumed to be abandoning the
extended-improvement tie-break — flagged as the round's quality risk (§7.1).

**The central finding — the plan's premise is false in the favourable
direction.** The tie-break is compatible with an O(n) screen. Maintaining, per
out-of-selection point, `min1[s]` (distance to nearest selected) and `arg1[s]`
(the vertex realising it), with strict-`<` updates, keeps `arg1` a true witness
under FP ties; then for a candidate drop, `arg1[s] != drop` means the post-drop
minimum IS `min1[s]` (the witness survives, and `min` re-selects the same
matrix cell), and only the ~(n−m)/m candidates whose witness is being dropped
rescan in O(m). Every candidate still gets its exact `nd` in the same ascending
order, so the tie-break fires on the same candidates with the same values and
the chosen swap — hence the whole trajectory — is bit-identical. What the
plan's §3 actually rules out is *sub-linear* candidate skipping, not O(1)
exact enumeration. That lands the round in the plan's own §1 fast lane
("a change that only reduces work per iteration can be verified as before"),
and its §7 decisions collapse: **tie-break kept exactly; `.Grasp_R` untouched
(parity test still the oracle); the manuscript's cached GRASP figure data
remains valid and no FurthestPoint re-pin is needed (§7.3 moot).**

**Harness (built first, per plan §5):** `drivers/grasp-frontier.R` (78291a3) —
mode (a) equal-plateau speed + objective identity at canonical n = 2000 shapes
the battery's n ≤ 500 grid never reaches; mode (b) equal `maxSeconds` quality,
plateau-off and canonical-512 sub-modes, 5 seeds/cell, mean AND worst per-seed
deltas, `iters`/`pr_calls` recorded so budget-bound phase-C skipping is
visible as mechanism. Real case: satellite (tc14) thinned to 2000 by
farthest-first, mirroring the canonical coreset pipeline.

**T-016 — LS nearest-selected screen (6ab0cb9):** dominant scan O(n·m) → O(n)
per pass; maintenance under a swap is one contiguous column read plus rescans
of the dropped vertex and its stranded witnesses. Verified: battery 1458/1458
(indices, objective, `iters`, `pr_calls`), R↔C++ parity 324/324, suite
4893/0/6. Clean-build timing (n = 2000, k = 100): plateau 8 0.47→0.22 s
(2.1×), 64 2.21→0.62 s (3.6×), 256 10.49→2.58 s (4.1×); k = 10 flat — the
T-015 small-m regression trap avoided. Instrumented counters (plateau 256):
screen rescans 1.6% of candidates (≈ 1/m prediction), evictions ~30/swap,
tie-arm fires 0.8% with 94% memo hits — all new branches exercised.

**T-016b — LS header fold (eed52c7):** `dstar` IS the running `gmin`; the pair
count becomes a reset-on-new-min counter inside the `di` sweep; the ≤2
`objective_of(sel \ witness)` rescores become one lazy shared exclusion pass.
Worth a small consistent win only where the header binds — plateau 256
2.58→2.47 s, every T-016b run beating every lever-A run — flat elsewhere;
kept on the T-004 precedent. **Method note:** a first single-chain best-of-3
showed a phantom +36% regression at plateau 8 that did not reproduce; this box
throws ±20–35% spikes on sub-second cells, so keep/reject timing calls were
decided by interleaved min-of-N across separate invocations (the skill's
"never a single run" guard, applied via interleaving).

**T-017 — PR-walk cross-step maintenance (a5d1f48):** the walk's candidate
lists, global-min-edge witness (per-member nearest-other summary `edi`/`earg`
under the same witness invariant), and per-add-candidate `NearTwo` (now
carrying an `arg2` witness so eviction is detectable) all persist across
steps; witness drops rescore via the shared lazy exclusion sweep. Per-step
unconditional O(m²) → O(m) plus rare rescans. Battery 1458/1458 incl.
`pr_calls`. Interleaved min-of-3: plateau 8 0.20→0.12 s (1.7×), 64 1.1×,
256 flat — the expected profile for the PR-bound low end.

**Round result (cumulative vs PR3 tip fde54d9, n = 2000 k = 100):** equal
plateau 8/64/256 = **3.9× / 4.3× / 4.4×** (timing minima; frontier sums say
3.8×/4.2×/4.5×), satellite coreset 2.6–3.2× across capture rounds, k = 10
1.0–1.2×, n = 500 k = 50 ~2×. **Objectives IDENTICAL on every equal-plateau
frontier cell.** Equal-budget: **119/120 seeded cells win or tie**, mean T_k
+0.14% to +1.64% at k = 100 with the mechanism visible (4–5× the iterations
per budget; path relinking reached inside budgets where the baseline got
none). The single loss (k = 10, canon, 0.25 s) did not reproduce: a targeted
3×5-seed rerun hit the best objective 15/15 on the new build while the
baseline itself dropped seeds under its own budget noise.

**Phase split after (instrumented, seed-42 instance, plateau 8/64/256):**
construct 0.015/0.056/0.351 s, LS 0.084/0.228/1.335 s, PR 0.064/0.045/0.044 s
(PR was 0.157/0.100/0.138 before T-017). **Next floor: the local search still
dominates both ladder ends (52%/69%/77%).** Its remaining mandatory term is
the per-pass m²/2 pair sweep that the extended-improvement pair count
requires (di + count + witness now share that one sweep), plus the O(n·m)
per-call summary init for PR-called LS (a construction `g`-handoff would
remove it for phase-B calls only). §8's construction-binds prediction is not
yet reached (construction 9–20%); its lazy-heap idea cannot beat
construction's own O(n)-per-step update pass and was not attempted.

**Findings placement (per the /profile skill):** all three levers were fixed
in-round → no GitHub issues, and no new rows in the frozen `findings.md`
archive; this entry and the commits are the record (T-016/T-017 are
log-narrative labels, not issue ids; issue filing is in any case unavailable
while `ms609-agent` is suspended).

**Timing policy:** all ratios are Windows-local *relative* figures (clean
builds, interleaved minima). No Linux/gcc cross-check this round — the levers
are op-count reductions, not CRT artifacts — but per the skill, canonical
wall-clock and any mission-wide claim await Hamilton.

**Docs:** NEWS 0.0.0.9007 states plainly that selections are *unchanged*;
DESCRIPTION bumped; baselines.md Area 3 refreshed to the round-4 build.

**Cleanup:** no VTune result dirs created; scratch libraries and both
instrumented builds live under the session scratchpad, outside the repo;
no instrumentation ships (verified by reverting `src/grasp.cpp` to the
committed state after each instrumented capture). last_focus unchanged
(targeted run continuing round 3, not a rotation slot).

**Coverage addendum (suite, not driver):** a cumulative branch-counter build
run against `tests/testthat/test-grasp.R` alone (4266 pass / 0 fail) shows
every branch added this round fires under the existing tests — LS screen
rescans 297,878, witness evictions 120,756, dropped-vertex fresh scans
26,968, exclusion pass 39,370 (both endpoint arms); PR exclusion pass 3,832
(both arms), `edi` fresh/evict/adopt 4,298/4,336/16,312, NearTwo witness
evictions 4,098 and O(1) inserts 4,276. No targeted test additions owed;
CI codecov is expected green on the existing suite.
