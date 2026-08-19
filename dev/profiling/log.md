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

---

## Round (area #6) — MaxMean RLTS tabu loop — 2026-06-16

Target: the per-iteration tabu inner loop in `src/maxmean.cpp` (user-named).
MaxMean is time-budgeted, so the metric is THROUGHPUT (iters/s) — more iterations
in the fixed budget = better solution quality.

Driver: `drivers/maxmean.R` (signed n=500 Type-I, 3 s, useRL=FALSE to isolate the
tabu loop). Baseline 786k iters/s.

profvis triage: skipped (the R wrapper is a coercion + symmetrize + one long
`.Call`; >99.9% of time is the pure-C++ kernel). Went straight to VTune.

VTune hotspots (`-g` build, result_maxmean): per-iteration self-time split ~69%
best-flip scan (lines 229–247) / ~31% O(n) P-array update (256–271). Hottest
lines were the two per-element divisions (233 `/(m+1)`, 237 `/(m-1)`, ~0.84 s) and
the P-update stores (260/268).

Micro-bench triage (the gate — VTune line shares are NOT a verified win):
- `bench-scan.cpp`: plain strength-reduction (precompute 1/(m±1), multiply) =
  **0.97×, REFUTED**. The scan is ILP/load-bound at n=500 (arrays L1-resident);
  the divide overlaps other work, so the "hot division line" was loads feeding it.
- `bench-scan2.cpp`: (a) monotonicity scan = **1.13×** (track max-p add / min-p
  remove with precomputed aspiration thresholds; divide only 2 survivors;
  best_flip+delta identical). (b) branchless P-update = **1.4×+** (zero the diag at
  the R boundary so the j!=u guard drops; bit-identical p-arrays).

Implemented both (T-012). End-to-end driver on the optimized -O2 install:
**786k → 1.16M iters/s (1.47×)** — better than the ~1.20× Amdahl estimate; solution
quality also improved (f 54.524→54.541, more iterations). Full `test_dir` green;
covr 100% on all new code (157/157 C++, 42/42 R) — the rewritten scan's new
branches are all exercised.

Correctness: branchless P-update is bit-identical (guard only skipped the
now-zeroed diagonal). Monotonicity scan changes only cross-type exact-δ tie-breaks
(measure-zero on continuous distances); all optima/consistency tests pass.

- status: area #6 OPTIMISED (T-012). The P-update is now branchless memory-stream
  (likely near AT-LIMIT — GCC -O2 doesn't auto-vectorize, but the loop is simple);
  the scan is monotonicity-reduced. A further round could test -O3/`#pragma GCC
  optimize("tree-vectorize")` on the P-update loop, or narrowing tabu_until to int
  for large-n bandwidth — both speculative, deferred.
- cleanup: result_maxmean_* and .vtune-lib-* deleted.

---

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

---

## Round — 2026-08-13 — Area 3: Grasp round 5 (to the ceiling, then across the cores)

**Question (user):** the project compares heuristics on a time/quality
frontier, so every implementation must be absolutely optimised — keep going
until the ceiling, parallelism authorised with core count controlled per the
`../TreeDist` convention. Design in ROUND5_PLAN.md; order was load-bearing
(the RNG change invalidates the battery baseline, so exact levers landed
first).

**T-018 — member summaries across passes (cde916a).** Round 4's declared
next floor — the per-pass m²/2 pair sweep — dies: each selected member
carries (mmin1, mmin2, witnesses for both slots, and mcnt = #partners AT
mmin1). dstar, its witness edge, the extended-improvement pair count
(halved mcnt sum over members at dstar — even and exact) and the post-drop
witness rescores (near_excl per member, O(m) per endpoint) all read off it;
a swap adjusts O(m) entries plus rare witness-eviction rescans, with a tied
departure decrementing mcnt BEFORE the entrant's insertion. All working
arrays moved into an LSScratch allocated once per Grasp_cpp call. Battery
1458/1458 bit-identical; interleaved min-of-3: plateau 8/64/256 =
0.12→0.09 / 0.50→0.37 / 2.25→1.68 s (~1.35× uniformly, every run beating
every round-4-tip run).

**Tie coverage hole closed (de87ef9).** Continuous random data cannot
produce two equal distances from different cells, so the tie arms
(duplicate minima, mcnt > 1, the tie-decrement, multi-pair dstar) were
unreachable by every existing test AND the battery. A lattice-points parity
test supplies exact ties in bulk; branch counters confirm the suite now
fires every T-018 branch (tie-decrement 9, duplicate-min insertions 681,
multi-pair dstar 204). Suite 4897/0/6.

**T-019 — construction scan fold (603c779).** The gmax/gmin/argmax scan
merges into the previous step's g-update pass (g[pick] = -Inf set first;
same values, same ascending order, same first-index ties); one standalone
scan seeds the loop; the dead final update is skipped; the RCL buffer is
reused. Bit-identical 1458/1458; ~1.09× overall (0.37→0.34, 1.69→1.55 s),
every T-019 run beating every T-018 run.

**T-020 — measured, dropped.** Post-T-018 the LS init (candidate summary +
member seed) is 9–11% of total wall (tLSinit 0.091 of 0.996 s at plateau
256). A construction-g handoff would save ≤ ~8% for real coupling; below
the bar. Serial split after T-019: LS passes ~64%, construction 22%,
PR 4%, LS init ~9% — AT-LIMIT-shaped (the screen's O(n) is mandated by the
tie-break's enumeration requirement, construction by its own sequential
update pass).

**T-021 — nCores-invariant parallel GRASP (6295fb8).** The round's one
deliberate trajectory change, landed last. Constructions consume uniforms
pre-drawn on the main thread (floor-index, bias O(size/2^53), fixed m draws
per construction; GRASP_BATCH = 32 is algorithm-defining since R-stream
consumption depends on it); workers (OpenMP, schedule(dynamic), per-thread
LSScratch pool) run construct + LS + objective with no R API; batches merge
in iteration order with all stopping rules at merge time; phase C
parallelises its pairs with a pair-order reduce. `mc.cores` (TreeDist
convention, default 1) → n_threads. **Contract verified: one seed,
nCores ∈ {1,2,8}, three shapes — bit-identical** (grasp-invariance.R).
R↔C++ parity holds over the whole grasp suite (4271/0 incl. the lattice)
because .Grasp_R mirrors the batch protocol draw for draw. Quality:
equal-plateau deltas −0.24%…+0.54% (mean +0.2%); equal-budget k = 100 up to
+1.6% mean T_k with PR reached inside budgets the old build exhausted; two
losing seeds only at the tightest k = 10 quarter-second budget (−0.33%
mean). Scaling (16-core box, grasp-scaling.R): plateau 256 =
0.83/0.43/0.23/0.14 s at 1/2/4/8 threads (**5.9× at 8**), objective
constant across the row; per-iteration serial cost unchanged (~2.1 ms at
plateau 256). New battery baseline captured (c5c-battery.rds pattern) —
pre-T-021 captures are obsolete as anchors. The phase-C scan test was
retuned (k=6, es=8, alpha=0.3, n=40): the old shape's scan contained no
PR-improving case under the new stream in 400 seeds (replay verified exact
— 0 ceiling violations — before retuning).

**Environment trap (memory: makevars-pkgcxxflags-clobber):** the user-level
`~/.R/Makevars.win` assigns bare `PKG_CXXFLAGS =`, which empties every
package's own src/Makevars compile flags — OpenMP linked but never
compiled, silently serial. Agent builds redirect R_MAKEVARS_USER to a
scratch copy (ccache/-j8/-O2 kept, clobber dropped); the package ships
canonical $(SHLIB_OPENMP_CXXFLAGS) in PKG_CXXFLAGS and PKG_LIBS. Reported
to the user — their own local OpenMP builds (TreeDist included) are
affected.

**Round-5 result.** Serial: 2.25 → 1.55 s at the canonical plateau-256 cell
(1.45×, bit-identical to round 4's trajectory), on top of round 4's 4.4×.
Parallel: ×5.9 at 8 threads with seed-determinism independent of core
count — the strongest reproducibility property the manuscript comparison
can hold. Manuscript figure caches: invalidated ONCE by T-021 (the round-4
"no re-pin needed" no longer stands); future re-runs are core-count-proof.

**Ceiling statement.** Serial residue is at its structural floor: the
extended-improvement rule requires enumerating every candidate (O(n) screen
per critical drop), construction is bound by its own O(n)-per-step
sequential update, the tie-arm memo is the last m²/2 term and fires ~2×
per pass. The remaining axis is core count; real scaling curves and
canonical wall-clock belong on Hamilton, with the Linux/gcc cross-check
owed before any mission-wide claim (/profile skill).

**Timing policy:** Windows-local relative figures throughout (interleaved
min-of-N per the round-4 method note). **Cleanup:** no VTune dirs; scratch
libs and instrumented builds in the session scratchpad; instrumentation
reverted after each capture. last_focus unchanged (targeted continuation).

**Round-5 closing addenda.** Full suite on the final build with the retuned
test: **4898 / 0 fail / 6 known-err**; invariance additionally confirmed at
nCores = 16 (the box's maximum). Coverage caveat: unlike T-018 (whose
branches were counter-proven under the suite), T-021's mid-batch budget-skip
and phase-C worker-skip branches need the wall clock to expire inside a
batch — inherently flaky to test deterministically — and the two floor-index
clamps are defensive-unreachable by design; codecov may flag those lines
(`# nocov` markers are a reasonable PR-time answer). The OpenMP build has
compiled and run only on Windows/MinGW so far — Linux/gcc is unexercised
until a push triggers CI, and real scaling curves belong on Hamilton.


## Round 6 — 2026-08-13 — Area 3 (Grasp): VTune certification + two exact levers

**Trigger:** user mandate to eke out every remaining ounce before the PR.
Rounds 4–5 optimised via timing + branch counters; the shipped 0.0.0.9008
kernel had never been VTuned. First hotspots collection (symboled -O2 -g
-fno-omit-frame-pointer build WITH $(SHLIB_OPENMP_CXXFLAGS) — the skill's
stock recipe would have stripped OpenMP and profiled a kernel without its
parallel regions) on the canonical serial shape (n=2000 dim=10 k=100
plateau=256, 6 reps, 5.05 s bare): grasp_local_search 32.8%,
grasp_construct 19.1%, **count_pairs_within 10.7% + count_to_set_le 6.6%**
— the tie-arm, which round 5's ceiling statement had written off as
structural, was the single biggest addressable residue.

**Lever 1 — tie-break count derived, not counted (exact).** npc =
#pairs(rem ∪ {s}) ≤ nd splits into rem×rem and s×rem halves whose floors
are both already known exactly: base_z IS rem's pairwise floor, cross IS
s's floor to rem, and nd = min(base_z, cross). So the rem half is
(nd == base_z) ? #pairs in rem at exactly base_z : 0 — and that count is
pair_count − mcnt[drop] when base_z == dstar (O(1) off the T-018
summaries), else counted from the rows of the ≥2 members achieving the new
floor (read off near_excl, O(m) each, once per drop, lazily). The s half is
(nd == cross) ? #v ∈ rem at exactly cross : 0, a direct O(m) row count run
only for win-or-tie candidates. Same integers from the same exact double
comparisons ⇒ trajectories unchanged.

**Lever 1a — REJECTED by measurement: maintained candidate tie count.** A
cnt1[s] (#selected at exactly min1[s]) maintained through init/swap would
make the s half O(1); it measured 1.17/1.11/1.19× on the k=100 ladder but
0.82× at k=10 over 40-rep interleaved runs — the maintenance streams a
second column (col_drop) per swap, a cost ∝ n that small-m tie-arms never
repay, and even at k=100 it underperformed the stateless variant
(1.25/1.21/1.32×). Dropped; the O(m) s-half count kept.

**Lever 2 — construction→LS g-handoff (T-020 revisited, exact).**
grasp_construct now keeps g in W.min1 with a pick-order witness in W.arg1,
folds the final pick's column back in with one lean pass (T-019 had skipped
it as dead), and the seeded local search skips its m-column candidate-init
re-read. min accumulates exactly, so values are bit-equal to the old
column-order init; the witness settles on first-pick-to-reach rather than
lowest-index under ties — a different but valid witness that routes the odd
rescan differently while every value read is identical. Phase-C local
searches (path-relinked starts, no construction) keep the plain init.

**Verification (each variant separately):** battery 1458/1458 bit-identical
vs the clean 0.0.0.9008 build + R↔C++ parity 324/324, both levers; full
suite on the final build 4898 / 0 fail / 6 known-err / 2 skip.

**Measured (interleaved min-of-3, single-threaded, this box):**
| shape | 0.0.0.9008 | round 6 | speedup |
|---|---|---|---|
| n=2000 k=100 plateau=8   | 0.14 s | 0.11 s | 1.27× |
| n=2000 k=100 plateau=64  | 0.20 s | 0.15 s | 1.33× |
| n=2000 k=100 plateau=256 | 0.84 s | 0.61 s | 1.38× |
| n=2000 k=10 plateau=64 (40-rep) | 1.10 s | 1.02 s | 1.09× |
| n=500 k=50 plateau=64    | 0.04 s | 0.04 s | ~1× |

**Post-round VTune (same driver, 7 reps):** grasp_local_search 44.2%,
grasp_construct 25.4%, grasp_path_relink 3.4%; the two counting functions
are gone from the profile. **Revised ceiling statement:** every remaining
hot loop is a single ordered pass the semantics require — the O(n)
candidate screen per critical drop (extended-improvement rule), the O(n)
per-step construction update (each step must stream the picked column), the
O(m²/2) member-summary seed per LS entry, and the O(n) col_add maintenance
per swap. No further sub-linear substitution is visible; serial is at its
structural floor with VTune evidence this time. Remaining axes: cores and
Hamilton.

**In-round fixes, no issues filed** (skill rule); driver added:
drivers/grasp-vtune6.R. Baselines refreshed (regress against round-6 rows).
Cleanup done: result_grasp6*/ and .vtune-lib-* deleted post-round.
Environment note: the user deleted the bare PKG_CXXFLAGS line from
~/.R/Makevars.win this round (memory updated) — package OpenMP flags now
flow to agent AND user builds without the R_MAKEVARS_USER workaround; the
symboled VTune build still needs one (adds -g -fno-omit-frame-pointer via
PKG_CXXFLAGS without losing $(SHLIB_OPENMP_CXXFLAGS)).
last_focus unchanged (targeted continuation of area 3).



## Round 7 — 2026-08-13 — Area 1 (FarFirst): argmax fold + nCores-invariant parallel kernels

**Trigger:** user mandate ("direct the same level of effort at ... the
FarthestFirst solver"); area files changed 2026-07-08 > last profiled
2026-06-10, so the rotation was due regardless. Prior ground re-read first:
round 1 priced the R glue at < 0.2 ms (T-003 refuted "fold into C++"),
T-001/T-002 shipped free-T_k + squared space, T-004's column reorder left
the points pass "AT-LIMIT" — but the argmax fold and parallelism were never
tried on these kernels. VTune (symboled OpenMP build, farfirst-vtune7.R:
matrix N=6000 n=3000; points dim {2,10}; N=60,000 case): EuclidColSqInto
1.63 s + MaximinFromPoints self 0.87 s (the separate argmax/pmin passes) of
~2.9 s kernel time; matrix kernel 0.23 s, bandwidth-bound.

**Lever 1 — fused greedy step (exact; the T-019 pattern).** Both kernels
seed (best, best_val) with one standalone scan, then ride each step's
argmax on the update pass — strict >, ascending i, masked −Inf entries
can't win, so the first-maximum tie rule is preserved read-for-read. The
points kernel additionally fuses the squared-column fill into the same
sweep: dimension 0 writes (no zeroing pass), middle dimensions accumulate
in the same j-ascending order (identical partial doubles), the last
dimension finishes in a register, merges, and tracks the max. Final step's
dead update skipped (nothing reads min_dist after the last pick). dim ∈
{0, 1} handled (0 = degenerate all-zero column, kept defined).
**Measured** (interleaved min-of-3): points dim=2 **1.43×** (N=6e3) /
**1.49×** (N=1e5); dim=10 1.11× / 1.14×; matrix 1.04× (bandwidth-bound, as
VTune said). Battery 925/925 bit-identical + cross-path identity intact.

**Lever 2 — mc.cores parallelism, results invariant to thread count.**
No RNG anywhere in these kernels, so exactness is chunking-trivial:
elementwise pmin/accumulation each computed by exactly one thread in the
serial order; argmax reduced per contiguous chunk then merged ascending
with strict > (= global first maximum); RowSums/RowSqSums keep each row's
long-double accumulation on one thread in j order; Diameter merges
column-chunk winners in ascending chunk order (column-major first-wins
preserved). Greedy pass engages at nPts >= 32768 (GONZ_PAR_MIN — below
that the per-step barrier outweighs a ~30 µs sweep; pure tuning constant,
results identical either way). Anchors parallelise at nPts > 64.
`.NThreads()` (utils.R) reads mc.cores once per wrapper; Grasp() refactored
onto the same helper. Oracle path stays serial (calls back into user R).
**Measured** (8 threads): N=1e5 pass 2.5× (dim 2) / 4.5× (dim 10);
O(N²) anchors 5.3× (compute-bound, near-linear as expected); N=6e3 cells
unchanged (below threshold, by design). Serial cost of the parallel build:
within noise (0.96–1.00×). **Verification:** battery bit-identical at
mc.cores 1 AND 2 (925/925 each, anchors' parallel paths exercised);
farfirst-invariance.R INVARIANT over nCores {1,2,8} × 5 cases including
two above-threshold passes and the seeded default ensemble; two
mc.cores-invariance tests added to test-farfirst.R (CRAN 2-core cap
respected).

**Round-7 result.** Cumulative on the round's cells: points N=1e5 dim=10
0.63 → 0.125 s (5.0×, 8T), dim=2 0.203 → 0.057 s (3.6×, 8T); anchors
ensemble N=6e3 0.73 → 0.14 s (5.2×, 8T); serial-only gains 1.04–1.49×.
Every gain is selection-preserving: 925-case battery + 1620-style
cross-path identity + suite. Area 1 stays AT-LIMIT serially (matrix path
bandwidth-bound; points pass now one fused sweep per step); remaining
axis is cores, which now scale. Windows-local relative figures
(interleaved min-of-N); Linux/gcc first exercised by this branch's CI.

**In-round fixes, no issues filed** (skill rule). Drivers added:
farfirst-battery.R (925-case old-vs-new + cross-path), farfirst-timing.R,
farfirst-invariance.R, farfirst-vtune7.R. Baselines refreshed (area 1
round-7 rows). Cleanup: result_ff7/ and .vtune-lib-* deleted post-round.
last_focus unchanged (user-targeted round on area 1).



## Round 8 — 2026-08-13 — Area 3 (Grasp): memory-layout stab, REFUTED — floor certified on the second axis

**Trigger:** user-requested "one last stab" after PR #2 merged. The branch
was first rebased onto post-#2 main (clean, no conflicts; interleaved
re-measurement matched the round-6 baseline rows: 100/150/590 ms on the
plateau 8/64/256 ladder — no regression from the rebase).

**Candidate lever — symmetric-transpose reads (bit-identity preserved).**
Round 6 certified the algorithmic passes but never the *orientation* of the
matrix reads. Fresh line-level VTune on the rebased tip put three stride-n
gather loops at ~29% of the 4.3 s kernel: the `arg1[s] == drop` cross
rescan (grasp.cpp:446, 0.76 s ≈ 17% alone), the post-swap candidate rescan
(:530, 0.40 s) and the floor-achievers count (:465, 0.11 s). `d` is
symmetric by documented contract, so each `d(s, v)` can be read as
`d(v, s)` — one contiguous walk down a single column instead of a gather
across the m selected columns, the same doubles bit for bit. Implemented at
all three sites (plus the s-half tie count sharing :446's shape);
verified bit-identical before timing: battery 1458/1458 vs the rebased-tip
build, R↔C++ parity 324/324, iters/objective identical on every ladder
shape.

**Measured — SLOWER at every shape; reverted.** Interleaved min-of-9
ladder + 3×40-rep k=10 cell, single-threaded:
| shape | base | transposed | ratio |
|---|---|---|---|
| n=2000 k=100 plateau=8   | 0.10 s | 0.14 s | 0.71× |
| n=2000 k=100 plateau=64  | 0.15 s | 0.20 s | 0.75× |
| n=2000 k=100 plateau=256 | 0.59 s | 0.86 s | 0.69× |
| n=2000 k=10 plateau=64 (40 calls, seeds 1:40) | 1.58 s | 3.33 s | 0.47× (worst rep 0.33×) |
| n=500 k=50 plateau=64    | 0.04 s | 0.04 s | ~1× |

**Why the profile misled:** the "gathers" all land inside the m *selected*
columns — an ~m·8n-byte working set (1.6 MB at the canonical shape) that
every rescan and every swap re-touches, so it stays cache-resident and the
VTune time on those lines is hot-set latency, not DRAM stalls. The
transposed orientation reads each rescanned candidate's *own* column —
~n/m fresh 16 KB columns per critical drop and per swap with no reuse
across candidates or swaps (~3 MB/swap of cold traffic at k=10, where the
regression peaks at 2–3×). The shipped orientation is the cache-optimal
one of the two; no hybrid is worth carrying (the transposed arm wins
nowhere measured).

**Verdict: AT-LIMIT on the memory axis too.** Round 6 certified every hot
loop as a semantics-mandated single pass (instruction axis); round 8
certifies the access orientation of those passes (memory axis) by refuting
the only alternative layout. Serial Grasp is at its floor with measured
evidence on both axes; remaining axes stay cores + Hamilton. No code
shipped (edit reverted), so no NEWS entry or version bump; no issues filed
(in-round refutation, skill rule). Post-rebase full suite green before
push. Drivers: committed grasp-timing.R/grasp-battery.R plus a scratch
40-rep k=10 loop (seeds 1:40 — note this draws a different workload mix
than round 6's fixed-seed 40-call row, so the two k=10 totals are not
comparable; ladder rows are). Cleanup done: result_grasp8/, symboled and
scratch libs deleted post-round.
last_focus unchanged (targeted continuation of area 3).

## Round 9 — 2026-08-13 — Area 1 (FarFirst): validation, block sweeps, anchor scans — serial floor certified

**Trigger:** user mandate ("optimize the performance of FarFirst to floor …
keep working until unimprovable"). Fresh line-level VTune on the round-7 tip
(farfirst-vtune9.R, which adds the anchor primitives round 7 never
line-profiled), plus a components triage of the matrix cell.

**The triage finding that set the round:** the matrix timing cell spent
220 ms/call of which the greedy kernel was **10 ms** — `.AsDistMatrix`'s
`anyNA(d) || any(!is.finite(d))` guard cost 200 ms (two full-size logical
intermediates, ~288 MB of allocation at N = 6000) and had never been triaged
because round 7 profiled only the kernels. Every shipped lever, in order,
each gated on the 925-case battery + cross-path dim sweep (88 cells) +
nCores-invariance BEFORE timing (interleaved min-of-3, serial):

- **A — AllFinite_cpp** (src/utils.cpp): single-pass allocation-free finite
  scan (branch-free exponent-mask OR-reduction, integer ops so it needs no
  FP reassociation licence; OpenMP reduction past 2^20 elements,
  order-independent by construction). `.AsDistMatrix` swaps the idiom for
  it, so DropAdd/Grasp/KCentre matrix intake gains too. Matrix cell
  230 → 47 ms (**4.8×**).
- **B — blocked dimension sweeps** (SweepBlock, maximin_points.cpp): the
  fused update processes up to four dimensions per sweep with the running
  squared sum in a register — dim ≤ 4 never touches the scratch column;
  larger dim cuts dc traffic from one RMW per middle dimension to one per
  block. Same left-associated per-element chain, store/load of a double is
  exact ⇒ bit-identical. Points cells **1.35–1.40×** across dim 2/10,
  N 6e3/1e5.
- **C — Diameter triangle + squared-space guard**: the column-major first
  max of a symmetric matrix always lies in the strict lower triangle
  (pair {a,b}, a<b: index b+aN < a+bN), and relative order among lower
  cells is preserved ⇒ triangle-only scan is tie-identical. The scan runs
  in squared space with `best == sqrt(bestSq)` as invariant; sqrt only on
  running-max candidates (a sq ≤ bestSq pair can never win the reference's
  strict >; a sq > bestSq pair gets the reference's own sqrt-space
  comparison, because two distinct squares can round to the same sqrt).
  Parallel chunks re-balanced by pair count.
- **D — RowSums/RowSqSums serial pair-halving**: lower-triangle sweep adds
  each pair's distance to both endpoint rows; every row still receives its
  contributions in ascending-index order, the mirrored distance is the same
  double ((-x)·(-x) == x·x), and the dropped diagonal contributed an exact
  +0. Parallel path unchanged (per-row full sweeps) — pair-halving would
  interleave rows' update order across threads. C+D: points anchors
  ensemble 730 → 390 ms (**1.87×**).
- **E — matrix anchor scans** (MatrixOffDiagMax_cpp, RowSqSumsFromMatrix_cpp,
  maximin.cpp): the R idioms copied the N × N matrix (diameter: dOff with
  −Inf diagonal) or materialised d² (rownorm). Full-matrix scans — NOT
  triangle: `.AsDistMatrix` admits asymmetric matrices, where an
  upper-triangle cell can win — verified against the R idioms on
  asymmetric and tie-dense inputs at 1/2/8 threads. Matrix anchors
  ensemble 620 → 200 ms (**3.1×**).
- **G — fused row-aggregate sweep** (RowSumsSqFromPoints_cpp): an ensemble
  whose anchors span both aggregate families (anti_medoid/medoid/rowsum
  need sums; rownorm needs squared sums) previously ran two full pair
  sweeps, each with its own 18M sqrts at N = 6000. One fused sweep fills
  both: each accumulator receives the identical summands in the identical
  order as its dedicated kernel (verified fused == dedicated at 5 sizes × 3
  thread counts, plus 108 both-family ensembles old-vs-new). Points
  anchors 390 → 290 ms (**1.34×**).

**Refuted in-round (no code shipped, no issues per skill rules):**
- **F — matrix-pass parallel threshold.** GONZ_PAR_MIN = 32768 is
  unreachable for the matrix kernel (8.6 GB matrix), so its parallel path
  is dead code in practice. Lowering it to 4096 in a scratch build and
  timing N = 16384 (2 GB): 8 threads **60 ms vs 30–50 ms serial** — the
  ~2000 per-step barriers dominate ~15 µs sweeps. Verdict: the matrix pass
  is serial at every RAM-feasible size; constant kept (aligned with the
  points kernel), comment + Rd updated to say so.
- **H — register-banded matrix row aggregates.** Hypothesis: the x87
  fld/fadd/fstp round-trip into the long-double accumulator array dominates
  RowSqSumsFromMatrix_cpp, so banding 4 rows' accumulators into registers
  should win ~1.5×. Measured: banded == accumulator-array == base-R
  rowSums == 60 ms at N = 6000 (and the cell A/B read 0.95–1.04×).
  Reverted; the simpler column-outer orientation ships. The matrix
  row-aggregate loops measure identical in every orientation tried — at
  their floor.

**Round-9 result (definitive interleaved A/B vs the round-7/8 tip, serial,
scores identical on every cell):** matrix N=6e3 k=3000 218 → 42 ms
(**5.1×**); points passes 1.27–1.38× (d2/d10, N 6e3/1e5); points anchors
730 → 280 ms (**2.6×**); matrix anchors 930 → 200 ms (**4.7×**, new
timing cell). 8 threads: N=1e5 pass d10 125 → 85 ms, d2 57 → 47 ms,
points anchors 140 → 60 ms, matrix cell 42 → 28 ms (validation scan
parallelises). Remaining serial components all measure at memory/semantics
floors: AllFinite and OffDiagMax stream at ~14.4 GB/s (single-core DRAM
ceiling), the pass's per-element argmax merge is a semantics-mandated
single pass (round 7), the only alternative access orientations measured
equal (H) or slower (F), and round 1's deferred Gram/gemv reformulation
stays excluded-by-contract (it rounds differently, failing the bit-identity
gate before any timing). Area 1 → AT-LIMIT serial, certified on the
instruction axis (rounds 1+7) and now the memory/validation axis (this
round); remaining axis: Hamilton wall-clock.

Verification stack: battery 925/925 bit-identical vs the round-7 tip;
cross-path identity at dims 1–11 × k {2,60,250,499} (88 cells, certifies
every SweepBlock instantiation); matrix-anchor helpers exact vs the R
idioms on asymmetric + tie-dense inputs; fused-sweep kernel == dedicated
kernels and 108 both-family ensembles old-vs-new; farfirst-invariance.R
INVARIANT over nCores {1,2,8}; full suite green (see commit). New tests:
AllFinite semantics + parallel branch, block-shape cross-path sweep
(test-farfirst.R), matrix-anchor scans + fused sweep (test-seed-score.R).
Drivers: farfirst-vtune9.R added; farfirst-timing.R gains the
ens-anchors-matrix cell. man/DropAdd.Rd re-synced with its roxygen source
in passing (the Concision commits shortened the source but left the Rd
stale). Cleanup: result_ff9 dirs and scratch libs deleted post-round.
last_focus unchanged (user-targeted round on area 1).

**Post-push addendum (same day) — block width swept, floor closed.** The
sweep width 4 was initially an un-swept tuning constant; a scratch-only
B = 6 variant ({6,4} split at dim = 10 — identical left-associated chain,
verified same score and indices) measured 390–400 ms vs 400 ms on the
N=1e5 dim=10 cell: within noise, refuted. Width 4 ships. Remaining
unpursued micro-levers, recorded so a future round need not re-derive
them (estimated ceilings all ≤ ~1.2× cell-level, below this box's
±20–35% noise floor for reliable single-cell verification): buffered
column-streamed fills for the serial pair sweeps (arithmetic would SIMD
but the scalar errno-guarded sqrt chain remains), banded pair
accumulators for the halved row aggregates (x87 RMW per acc[j] survives
banding — see lever H's null result), and a fused matrix
RowSums+RowSqSums pass (saves one 288 MB stream, ~20 ms of a 200 ms
cell).

## Round 10 — 2026-08-13 — Area 2 (DropAdd): blocked fills, fused-pass mc.cores; matrix kernel certified serial

**Trigger:** user mandate (continue "…then of DropAdd… until both are
unimprovable"), and PR #2's merge (oracle path) plus 2026-07/08 kernel
changes made the area stale (last profiled 2026-06-10). Branch
perf/dropadd, stacked on perf/farfirst for `.NThreads()` + the OpenMP
Makevars. New trajectory-identity battery (dropadd-battery.R, 295 cases):
kernels called directly with `want_trace = TRUE`, so every case pins the
full drop/add sequence, objective, secondary and iters — over Euclidean,
TIE-DENSE integer-valued, and asymmetric matrices, both m regimes, default
and explicit seeds, dims 1/2/7/10 — plus two within-build invariants
(matrix ≡ points trajectory at a shared explicit seed; the pure-R oracle
twin ≡ the matrix kernel), kept green throughout.

**VTune (symboled, kernel-direct):** points kernel EuclidCol self-time
~1.6 s of 2.2 s (per-point strided gathers + sqrt at all four column-fill
sites); matrix kernel's top cost the need_recompute branch (~37% at
m = n/2, scattered K×m reads), then the ADD argmax ~10%; the T-005a seed
row-sums stream at DRAM bandwidth (certified at floor). The wrapper's
`.AsDistMatrix` guard resurfaced as the matrix intake cost — already fixed
by round 9's AllFinite lever on the parent branch; attributed, not
re-fixed (timing cells call the kernels directly).

**Shipped levers (each gated on the 295-trajectory battery before timing):**
- **P2 — blocked fills + sqrt at consumption** (dropadd_mf.cpp): the four
  whole-column fills become contiguous dimension-blocked squared sweeps
  (round 9's SweepBlock shape; identical left-associated chains), and each
  element's true distance is sqrt'd as the record loops consume it — the
  fills stay pure SIMD streams, no separate sqrt pass re-streams the
  column, and every distance is sqrt() of the identical squared double.
  (First cut used a separate sqrt pass: 1.09–1.12× at dim 10 but **0.93×
  at dim 2** — the extra column round-trip; the consumption form fixed the
  regression.) Points cells **1.08× (d2) / 1.18× (d10 search) / 1.26–1.27×
  (d10, m = 10⁴)**; matrix kernel untouched.
- **mc.cores for the points kernel** (fused pass regions): each
  per-iteration pass — column fill + sqrt + record updates, with the ADD
  argmax riding the drop pass and the construction update pass — splits
  into one contiguous chunk per thread inside a single parallel region
  (two barriers per iteration). Exactness: every element is computed by
  one thread with the identical chain; chunk argmax winners merge in
  ascending chunk order (an earlier chunk's winner has a smaller index
  than everything after it — first-max preserved); count-zero elements are
  skipped in the pass and merged after the recompute with the explicit
  (min_dist, sum_dist, smaller-index) triple — a lexicographic maximum,
  order-independent; chunk-local need_recompute lists concatenate in
  ascending chunk order (= the serial push order); the objective loop and
  self-record gathers stay serial (the secondary's long-double sum runs in
  FIFO buffer order — order-load-bearing). Engages at n ≥ 16384
  (DA_PAR_MIN). **Measured at 8T:** construction m = 10⁴ **2.3×**
  (1400 → 610 ms), search 1.3–1.5× kernel-direct, **1.75×** end-to-end
  through the public API at n = 2e4 m = 2000; d2 positive (1.2×).
  nCores-invariant {1,2,8} through DropAdd() on three above-threshold
  cells + a 2-core CRAN-capped test; above-threshold full-trajectory
  identity (drops/adds sequences) verified 1T ≡ 8T on five cases
  including matrix n = 16384.

**Refuted in-round (no code shipped, log-only per skill rules):**
- **M1 — matrix fused argmax** (search + construction): bit-identical
  (295/295) but **flat** — n=4000 cells 0.90–1.05×, n=16384 construct
  1.11×, search 1.00×. The record arrays are cache-resident at every
  feasible matrix size, so the fold trades instructions for nothing
  (FarFirst's fold won because its arrays streamed from DRAM). Reverted.
  The same logic ships inside the points kernel's parallel path, where it
  is what makes the drop pass a single region — serially it measures flat
  there too (1400 vs 1460 ms construct, within noise the rest).
- **Matrix-kernel threading**: fused-region chunking at n = 16384
  (2 GB matrix, both m regimes) measured **flat at 8T** (230/550 ms
  unchanged) — each pass streams a fresh matrix column from DRAM, which
  one core already saturates. With M1 flat and threading flat, the matrix
  kernel is certified serial-at-limit: its per-iteration cost is the
  mandatory column stream + cache-resident record logic.
- **Fills-only threading** (first experiment): 1.1–1.3× and **0.79× at
  d2** — post-P2 the record loops dominate the iteration, so
  parallelising fills alone was Amdahl-capped; superseded by the fused
  full treatment above.

**Oracle path (PR #2's R harness): AT-LIMIT by design.** Rprof with a
realistic vectorised colFn (n = 4000, dim 5): the user's colFn is ~93% of
sampled time; the harness (.DropAddPick / .DropAddApplyAdd / record logic)
is single-digit %. No optimisation warranted; the R twin stays frozen as
the semantic reference.

**Recorded leads (not pursued; anti-dup memory):**
- Lazy second-minimum record for the matrix recompute branch (the top
  matrix cost, ~37% at m = n/2): (min2, count2) maintenance would replace
  most K×m rescans with O(1) pops; exactness surface is large (count
  semantics must reproduce the rescanned values exactly). Estimated ≤1.3×
  cell-level on the m = n/2 search.
- Symmetric-transpose reads for the matrix recompute/self-records (read
  d(x, s) from column x instead of column s): blocked by a documentation
  tension — DropAdd's Rd requires a symmetric matrix, but `.AsDistMatrix`
  documents that asymmetric matrices are silently accepted with d_ij and
  d_ji treated as independent, and this round's battery pins asymmetric
  trajectories. Resolving the contract either way is a maintainer call.
- K-row coordinate pre-gather for the points recompute branch (contiguous
  scratch for the K rows' coordinates): exact by copy; matters only at
  large m (recompute ≈ 8% of the m = 10⁴ profile).
- The parallel recompute merge argument — (min, count) partials merge as
  min-then-sum — is recorded here should the branch ever dominate.

**Baseline provenance:** the battery/timing baseline build is main at
0215450 (the post-#2 merge, perf/dropadd's original branch point);
regenerate with `Rscript dev/profiling/drivers/dropadd-battery.R
<out.rds>` against a scratch install of that commit, then compare any
later build with the two-argument form.

**Round-10 result.** Serial: points cells 1.08–1.27× (battery-exact);
matrix kernel unchanged and certified at its floor (seed row-sums at DRAM
bandwidth; fused argmax and threading both measured flat; recompute lead
recorded). Parallel: points path 2.3× construction / up to 1.75×
end-to-end at 8T, nCores-invariant, engaging at n ≥ 16384 — honest
DRAM-bound scaling, not FarFirst-class, and said so in the Rd. Area 2 →
OPTIMISED serial + mc.cores points path; matrix kernel AT-LIMIT. Drivers:
dropadd-battery.R, dropadd-timing.R, dropadd-vtune10.R,
dropadd-invariance.R. Cleanup: result_da10/ and scratch libs deleted
post-round.

## Round 11 — 2026-08-14 — Area 4: ExactMaxMin probe reduction (user-directed)

**User steer:** any provably-optimal witness is acceptable — only the proven
objective and the proof itself are load-bearing. That admits
verdict-preserving reductions (not alpha- or witness-preserving ones), which
round 5's bit-identity gate forbade.

**Baseline:** perf/dropadd tip fb93c2d; R-devel 4.7.0; highs 1.14.0.2 in a
scratch lib (highs had dropped out of the local libraries in an R upgrade,
so the exact tests had been skipping silently on this box). Grid = 14
ground-truth cases × k∈{2,4,6,10}: 56/56 proven+valid, 12.7 s single pass,
10.40 s interleaved min-of-5.

**Triage (Rprof, warmed):** highs solver_run 81–95% of wall on every
profiled cell (penguins k4 95%, ionosphere k10 93%, zoo k4 81%); Grasp warm
start 1.5–13.5%; R-side probe assembly noise. The process's FIRST IP pays
~0.5–0.6 s one-time highs init — it moves between cells across builds and
must be discounted in per-cell comparisons.

**Reduction-potential audit** (scratch audit-core.R): at the certifying
threshold — the smallest candidate above the heuristic warm-start value —
peeling the complement graph H (pairs ≥ λ) to its (k−1)-core plus a greedy
colouring settles 21 of the 60 real k∈{4,6,10} probes outright (k=2 cells
never probe: the warm start attains the diameter); where an IP survives,
greedy χ sits at k+0..3, so colour-class rows collapse the LP bound from
~n/2 to ~χ. tc21_spam (n=4601) is verdict-free at every k (cores 0/0/26);
tc14_satellite (n=6435) at k=4. An R-vector peel prototype cost 13–64 s at
n=6435 → the kernel had to be C++.

**Levers shipped (verified together as the round's delta):**
1. `ThresholdReduce_cpp` (src/exact_reduce.cpp): CSR adjacency, bucket
   (k−1)-core peel, component split, Welsh–Powell colouring, greedy-clique
   feasibility shortcut, and emission of the reduced per-component model
   (colour-class rows + cross-colour G-edges). Tie-breaks by vertex index
   throughout → deterministic output.
2. `.MaxISVerdict` reduce-first rewrite: empty reduction ⇒ infeasible with
   no IP; greedy k-clique ⇒ feasible with no IP; otherwise one IP per
   surviving component — colour-class clique rows (sum ≤ 1 each) plus
   cross-colour edge rows, same integer feasible set as all-pairs.
   Witnesses still validated against `d` independently of solver status.
3. Cap row sum(x) ≤ k on every component IP: the probe only asks α ≥ k, so
   feasible components stop at the first k-incumbent and infeasible ones
   prove bound < k without closing to their true α.
4. Setup-path vectorisation (post-IP-collapse Rprof on spam n=4601 showed
   the solve wall was ~50% `sort(unique(ud))` and ~30%
   `which(upper.tri)`+`arrayInd`, the IP entirely gone): direct
   column-major triangle index construction (`sequence`/`rep.int`) and
   radix sort + adjacent dedupe for the candidate grid — the identical
   `ui`/`uj`/`ud`/`cand` vectors. Stretch cells 1.26–1.38× (spam k4
   2.4→1.9 s, k10 2.6→2.0 s, satellite k4 5.4→3.9 s, interleaved
   min-of-3); grid cells unaffected (setup is µs at n ≤ 351).

**Formulation A/B (interleaved min-of-5, whole grid):**
- reduce + class rows, no cap: 5.88 s (1.77×) — but gaussians k6 0.66×,
  uniform k10 0.70× vs base.
- reduce + class rows + cap (SHIPPED): 5.42 s (**1.92×**); cap fixes
  gaussians k6 (1010→500 ms).
- reduce + all-edge rows + cap (classes deleted): 6.32 s — REFUTED: the
  colour-class rows beat highs' own presolve clique detection on nearly
  every IP cell (penguins 1.29–1.60×, density 1.55–1.58× relative);
  all-edges wins only uniform k10 (0.81×).
- Known trade: uniform k10 140→270 ms is the one cell the shipped
  formulation regresses (χ=11 vs k=10; highs' internal cover happened to
  beat Welsh–Powell there). Not chased: per-cell formulation switching
  would overfit this box.

**Verification stack:**
- test-exact-reduce.R (726 assertions, no highs needed for the reduction
  verdicts): kernel vs brute-force clique enumeration over random probes —
  verdict agreement, proper colouring, every G-edge covered by a class row
  or an emitted edge row, emitted rows exactly the cross-colour G-edges;
  targeted structures (peel cascade to empty, K222 colouring kill at the
  certifying-probe shape, K5 greedy stop at k, single-edge k=2,
  determinism); .MaxISVerdict IP-branch verdicts (C5 refute,
  bridged-triangle feasible where greedy misses, two-C5 refute-all,
  refute-then-feasible component order).
- Grid battery (drivers/exact-grid.R, NEW): 56/56 scores bit-identical to
  base, proven flags identical, witnesses independently valid. Scores stay
  bit-identical under witness changes because any optimal witness's min
  pairwise distance is the same double — cand[bestIdx] realised in `d`.
- Full suite on the shipping build; RNG contract test #4 (same-seed
  reproducibility) green — the reduction draws no random numbers.

**Result:** grid 10.40 → 5.42 s interleaved min-of-5 (**1.92×**); cell
extremes ionosphere k4/6/10 ≈5–7×, glass k4/6/10 ≈16×, five cells to ~0 ms
(verdict-free). **Reach:** spam n=4601 proven optimal at every
k∈{4,6,10} in ~2 s (1.9–2.0 s after lever 4) — the baseline build needed
201 s and 8.9 GB peak at k=4 (~100×, and the per-probe memory wall is
gone); satellite n=6435 k=4 proven in 3.9 s (baseline model ≈16 GB: not
attempted). Both stretch cases decide every probe combinatorially, so they
are formulation-invariant. The 683–990 k=10 targets on the shipping build:
breastcancer 4.0 s, pima 3.6 s, vehicle 20.6 s, vowel 172 s — vowel is the
one grid-family instance whose certifying IP stays genuinely hard
(core = n, χ = 13 vs k = 10).

**Recorded leads (unmeasured, not pursued):** per-probe O(n²) threshold
rescans (`ud < lambda` per probe); warm-start restart-count scaling at
large n (RNG-stream change now permitted); DSATUR in place of
Welsh–Powell (χ−k gap is 0–3 where IPs survive; a 1–2-colour improvement
would convert more probes to verdict-free); MIP start via highs_solve's
`start=` (wrapper supports it; incumbent value ≤ k−1 at the certifying
probe).

Status: Area 4 → OPTIMISED (rounds 5+11). Drivers: exact-grid.R (battery +
timing cells), exact-large.R (stretch, refreshed to the (k, d) API).
Cleanup: scratch libs lib-r11* deleted post-round; audit script + rds in
session scratchpad only.

## Round 10 — 2026-08-14 — Area 1 (FarFirst): restart-level parallelism shipped

**Trigger:** user question — a `FarFirst()` call with several restarts could
run one restart per core; would that pay? Rounds 7 and 9 certified the serial
floor and parallelised the *within-pass* sweep and the O(N²) anchors, but the
restarts themselves had never been an axis: `.ResolveEnsemble` ran them one
after another, and below GONZ_PAR_MIN each pass is serial, so a multi-start
call used exactly one core for the bulk of its work.

**Feasibility measured before implementing.** A scratch kernel (parallel-for
over seeds, per-thread `min_dist`, shared read-only input) calling the shipped
inner pass code, A/B'd against a serial loop in the same translation unit at
the same flags. Identity gated first (8 seeds × N=6000 × k=300, both paths,
nthreads 1/2/8, both arms == the shipped single-seed kernels). Scaling at
nSeeds=3: matrix 2.80×, points d2 2.95×, points d10 2.77× at 4 threads —
ceiling `min(nSeeds, cores)`, so ~90-98% efficient. The 2T rows land at
~1.4-1.5× because three equal-cost tasks on two threads take two rounds
(ideal 1.5×): intrinsic imbalance, not overhead.

**The bandwidth hypothesis was wrong, and the measurement said so.** The
matrix path was expected to saturate DRAM (three concurrent passes tripling
demand against the ~14.4 GB/s single-core ceiling round 9 measured) and cap
near 2×. It scales exactly like the compute-bound points path: a pass streams
only one 48 KB column per step, and saturation first appears at nSeeds=8
(5.08× not 8×), past the default. No lever was tuned to the wrong model
because the prototype ran before the implementation.

**Dispatch rule (measured, not assumed).** Above GONZ_PAR_MIN the existing
axis still wins — points N=1e5, k=1000, nSeeds=3: serial loop with 8T
in-kernel 365 ms (3.99×) vs restart-parallel on 3 threads 550 ms (2.65×),
because nS is the smaller number there. So the kernels take the
restart-parallel arm only when `seedThr > 1 && nPts < GONZ_PAR_MIN`, and
threads are clamped to the seed count (a pass cannot be split further;
un-clamped, 16 threads over 3 seeds *regressed* the prototype's points d2 cell
2.95× → 2.11× on idle-thread barrier joins). No nested regions: a hybrid
(nS seeds × threads/nS each) is the only thing that could beat 3.99× at large
N, and is unmeasured.

**Shipped.** `MaximinMultiFrom_cpp` / `MaximinMultiFromPoints_cpp`: the pass
bodies were extracted verbatim into `MaximinPass` / `MaximinPointsPass`
(caller-owned buffers, no R API), so the single-seed kernels and the
multi-seed kernels are the same code and cannot drift. All R allocation
precedes the parallel region. `.ResolveEnsemble` now de-duplicates the specs'
seeds, dispatches one batch, and maps back — replacing the per-driver
`gonzCache` environments, which existed for exactly that de-duplication.

**Round-10 result** (interleaved min-of-5, 8 threads, scores identical on
every cell): ens-default matrix nSeeds=3 80 → 40 ms (**2.00×**), points d10
207 → 72 ms (**2.88×**); nSeeds=8 174 → 50 ms (**3.48×**) and 547 → 123 ms
(**4.43×**). Serial unchanged (117→113 / 205→205 / 206→206 / 553→540 ms).
Anchor ensembles unmoved by design (matrix 106→106 ms): their O(N²) seeds
dominate and already parallelise, and their greedy passes are ~6 ms.

**Ceiling on the prize, recorded so it is not re-chased.** `nSeeds` defaults
to 3 and the documented quality knee is n ≈ 3-4, so a default call saturates
at ~3× and four cores exhaust it. Raising `nSeeds` buys wall-clock scaling but
little solution quality — `DropAdd()` remains the route to better solutions.
Area 1's remaining axis is Hamilton wall-clock.

Verification stack: battery 925/925 bit-identical vs the round-9 tip at
mc.cores 1 AND 2 AND 8; cross-path identity OK in every capture;
farfirst-invariance.R INVARIANT over nCores {1,2,8}; full suite green
(5146 pass, 0 fail). New tests (test-farfirst.R): multi-seed kernels ==
single-seed kernels at 1/2 threads over k ∈ {1,2,25}; their error branches;
the above-threshold serial-fallback arm at N=33000; default-ensemble
mc.cores-invariance (CRAN 2-core cap respected); and a seed-collision
ensemble certifying that de-duplicated dispatch still yields one
strategy_results record per label. Drivers: farfirst-timing.R gains three
ens-default cells. Baselines refreshed (area 1 round-10 rows).
last_focus unchanged (user-targeted round on area 1).

## Round 12 — 2026-08-14 — Area 4 (ExactMaxMin): the IP replaced by a clique
search; parallelism measured and declined

**Trigger:** user question — is there gain in parallelising ExactMaxMin, and
are further gains available there or in DropAdd? Baseline perf/exact tip
fa482b3 (round 11 shipped), R-devel 4.7.0.

**The parallelism question, asked of the thing that actually cost.** A phase
breakdown of the shipped solver put **one** HiGHS solve at 94-99% of wall on
every hard cell (ionosphere k10 0.72 s of 0.77 s; vehicle k10 17.75 s of
17.91 s) — every other stage was noise, so component- or probe-level
parallelism had nothing to divide. Re-posing vehicle's certifying model at
threads 1/2/4/8 and with `parallel = "on"`: **17.6-20.3 s, flat, no trend.**
HiGHS's MIP does not scale here, so the answer to the question as asked is no
— and that made the model itself the target.

**What the IP was being asked to do.** vehicle k10's certifying probe is
"prove no 10-clique exists in an 846-vertex graph of 31% density". A prototype
Tomita-style bitset branch-and-bound with a greedy-colouring bound answered it
in 4121 nodes; HiGHS needed 17.75 s to reach the same verdict (alpha = 9).
Across seven dumped certifying probes the prototype settled every one in
milliseconds, so the packing IP was removed rather than tuned.

**Levers shipped (verified together):**
1. `ThresholdDecide_cpp` (src/exact_reduce.cpp) replaces
   `ThresholdReduce_cpp` and the per-component IP: CSR build, (k-1)-core peel,
   component split, then a bitset depth-first search per component under
   Tomita's colour-sort bound, with the candidates of a node held as 64-bit
   words so a neighbourhood intersection is a word-AND pass. Component
   vertices are relabelled in descending surviving degree, so each node's
   greedy colouring follows Welsh-Powell order. Exhaustive: no clique found is
   a proof that none exists. Deadline shared across components; interrupt and
   clock checked every 1024 nodes.
2. Candidate grid restricted to the warm-start tail. Nothing at or below the
   heuristic's realised value is ever probed, so `cand` is built from the
   distances >= it and the search starts at `cand[1]`. The full sort was 53%
   of the wall at n >= 4601.
3. `TriangleAtLeast_cpp` / `EdgesAtLeast_cpp` replace `d[upper.tri(d)]` and
   the per-probe `ud < lambda` gathers. `ui`/`uj`/`ud` are gone: at n = 6435
   they held ~500 MB whose only consumer was the probe rescan, which is now a
   column-major scan emitting the H edge list directly.

**Root-branch parallelism prototyped, measured, NOT shipped.** Branch i of the
root loop needs candidates `{order[0..i-1]} & N(order[i])`, computable without
the sequential prefix removal, so root branches are independent. Measured on
satellite k10's five costliest probes at 1/2/4/8 threads:
- **Infeasibility proofs scale**: the 0.47 s probe → 0.14 s at 8T (**3.36×**),
  node count identical (8314) at every thread count — an exhaustive proof
  visits the same tree however it is divided, so the verdict stays
  deterministic.
- **Feasible probes do not**: 0.75-1.38×, and nodes inflate with threads
  (2259 → 6727 at 8T) as threads speculate on branches the serial descent
  never reaches. The serial search already finds a witness quickly.
- End to end this is worth ~1.2× on the one cell where search still dominates
  (satellite k10: decide 6.65 s of 14.8 s wall, spread over ~35 probes, none
  above 0.67 s), and it would make the *selection* thread-dependent — the one
  package-wide invariant this solver's relaxed contract does not already
  waive. Declined on that trade, not on the scaling.

**Verification:**
- Grid battery (56 cells): **0 mismatches** vs the fa482b3 baseline — score
  bit-identical, `proven` identical, every witness independently valid.
- Stretch cells: scores identical to baseline on all seven re-measured cells.
- test-exact-reduce.R rewritten for the new kernel: verdict vs brute-force
  clique enumeration over 240 random probes (exact agreement, not just
  verdict-preservation, since the search is exhaustive); the targeted
  structures; a **packing-IP oracle** at n = 30 via `highs`, past brute
  force's reach; and the triangle scans against the R idioms they replace.
- Full suite 6183 pass / 0 fail / 2 skip (Geo loader only).
- `highs` and `Matrix` namespaces confirmed unloaded after an ExactMaxMin
  solve; both stay in Suggests for ExactKCentre, MaxSum and the test oracle.

**Result** (interleaved min-of-3 for the grid; single runs for the stretch
cells, which are seconds to minutes):

| cell | base s | new s | speedup |
|------|-------:|------:|--------:|
| grid, 56 cells | 6.5 | 0.9 | 7.2× |
| breastcancer k10 (n=683) | 4.0 | 0.1 | ~40× |
| pima k10 (n=768) | 3.6 | 0.1 | ~36× |
| vehicle k10 (n=846) | 20.4 | 0.1 | ~200× |
| vowel k10 (n=990) | 171.4 | 0.2 | ~860× |
| spam k4 (n=4601) | 2.3 | 0.3 | 7.7× |
| spam k10 | 2.4 | 0.5 | 4.8× |
| satellite k4 (n=6435) | 4.5 | 1.0 | 4.5× |
| satellite k10 | — | 13.1 | new reach |

vowel k10 was round 11's one genuinely hard grid-family instance; it is no
longer hard. Round 11's recorded leads on the IP formulation (DSATUR, MIP
start via `start=`, per-probe threshold rescans) are **closed by deletion** —
there is no IP and no rescan left to improve.

**DropAdd, asked at the same time.** Round 10 (2026-08-13) certified the
matrix kernel serial-at-limit and shipped the points path's blocked fills and
mc.cores threading; nothing in that area has changed since, so the answer
comes from that record plus one new measurement prompted by PR #6's benchmark
comment. Interleaved local A/B of the four CI cells, **PR #6 head vs its
actual base `perf/farfirst`** (5 reps each, mc.cores unset as on the runner):

| cell | base ms | PR ms | change |
|------|--------:|------:|-------:|
| `FarFirst(20, d2000)` | 6.468 | 6.428 | +0.6% |
| `FarFirst(20, d2000, diameter)` | 5.752 | 5.842 | −1.6% |
| `DropAdd(20, d500, plateau=2000)` | 7.654 | 7.651 | +0.0% |
| `DropAdd(20, pts4000, plateau=1000)` | 145.9 | 125.3 | +14.2% |

The FarFirst cells cannot move: the only FarFirst-path difference between
those two trees is three deleted `// nocov` comment lines in `src/maximin.cpp`.
The runner's −12.81% and +12.88% on them are noise, as is the +9.09% on the
matrix DropAdd cell — round 10 shipped nothing for that kernel. The one real
movement is the points cell, which the runner understated (+2.73% vs +14.2%).
**The suite has no cell in the regime round 10 optimised**: mc.cores threading
engages at n >= 16384 and the blocked fills pay most at higher dim and large
m, while the cells are n = 4000, dim 8, m = 20. Remaining DropAdd leads are
unchanged from round 10: lazy second-minimum record (<=1.3×, one regime, large
exactness surface), K-row coordinate pre-gather (~8% of one branch), and
symmetric-transpose reads — **blocked on a maintainer decision** about whether
`.AsDistMatrix`'s documented acceptance of asymmetric matrices is a contract.

**Recorded leads (unmeasured):** warm-start quality at large n — satellite k10
spends ~35 probes because the heuristic bound sits far below the optimum, so a
stronger seed would cut both the gallop and the bisection, and it is now 8% of
that cell's wall in its own right; parallelising the O(n²) edge scan (21% of
satellite k10) is deterministic by construction if chunk fills are ordered.

Status: Area 4 → OPTIMISED (rounds 5, 11, 12). Drivers: exact-grid.R (highs
dependency dropped), exact-large.R. Cleanup: scratch libs and dumped probe
graphs deleted post-round.

## Round — 2026-08-14 — Area 6: MaxMean maxIter guard  [user: /profile MaxMean]

**Question:** the new `maxIter` feature (this branch) added a per-iteration
`iter < iter_cap` guard to the hot tabu inner-loop condition in
`src/maxmean.cpp`. Does that guard regress throughput (iters/s) — the metric that
sets MaxMean solution quality, since the search is time-budgeted?

**profvis / VTune skipped — deliberately.** The area-#6 round (2026-06-16) already
mapped this loop: >99.9% pure C++ from the Rcpp boundary, self-time ~69% best-flip
scan / ~31% P-array update, and a one-comparison edit does not move the hotspot
map. The question here is a *verified delta*, not a hotspot location, so it goes
straight to a controlled A/B (skill Step 6).

**Method — on-machine controlled A/B, both `-O2 -g` (build-symboled-lib.ps1):**
- A = HEAD (guard present: `while (depth < alpha_depth && iter < iter_cap)`).
- B = guard removed (`while (depth < alpha_depth)`), one-line edit, `src/` restored
  after the build. DLLs confirmed to differ; only `maxmean.o` changed, same flags.
- Driver `drivers/maxmean.R` (signed n=500, useRL=FALSE, 3 s), 6 interleaved reps,
  iters-in-3 s as the throughput proxy (interleaved min-of-N per this box's timing
  noise).

**Result — no regression.** A (with guard) measured *equal-or-faster* than B in
5/6 reps (median iters ratio A/B = 1.19; range 0.97–1.22). Since B does strictly
less work in the loop condition, "B slower" is logically impossible as a real
effect → the delta is codegen/layout, and the guard's true cost (one
predicted-taken compare per O(n) iteration, < 1 %) sits below the combined
run + layout noise. **The maxIter guard is throughput-neutral.**

**Absolute-throughput note (NOT a regression).** A clean plain `-O2` build
(default Makevars, no `-g`/frame-pointer) re-measured a stable ~0.97–0.98M iters/s
— below the T-012 round's 1.16M baseline but far above the 786k pre-optimisation
baseline, with the same solution output (f = 54.5405, |S| = 138; matches T-012's
54.541). The gap is cross-session machine variance (this box: ±20–35 % on timing
cells), NOT the maxIter change (the on-machine A/B isolates that) and NOT lost
optimisation (identical solution + >786k throughput confirm T-012 intact).
`baselines.md` left at 1.16M — an environmentally depressed reading is not a clean
baseline to overwrite it with. Aside: `-g -fno-omit-frame-pointer` itself costs
~9 % throughput (0.98M → 0.89M), noted for future symboled rounds.

**Verdict:** area #6 remains OPTIMISED (near AT-LIMIT). No issue filed — nothing
to fix; the maxIter guard is confirmed free.

- cleanup: all `dev/profiling/.vtune-lib-*` builds deleted.

## Round 13 — 2026-08-14 — Area 2 (DropAdd): symmetric-read recompute, on a
contract the user settled

**Trigger:** user decision. Round 10 recorded the symmetric-transpose lever as
blocked on whether `.AsDistMatrix`'s documented acceptance of asymmetric
matrices was a contract; asked, the answer was that it is not.

**The lever.** Each recomputed point `xx` needs `d(xx, S[j])` for every
surviving `j`, and a symmetric `d` offers both layouts. Which is cheaper is a
cache-line count: reading a row of each of the m columns of S touches m lines
(one useful double per line), while sweeping the m rows of column `xx` touches
`min(m, n/8)` lines. Column-major reads therefore win exactly when `m > n/8` —
and the measurement matched the model on both sides of it, so the kernel picks
per call. The construction's own-record scan transposes unconditionally: it
reads a single column either way, so the column form is never worse.

**Measured** (interleaved min-of-3, kernel-direct, n = 4000, objectives
identical on every cell):

| cell | base s | new s | |
|------|-------:|------:|---|
| m=10 construct | 0.0100 | 0.0100 | flat |
| m=10 search1500 | 0.0500 | 0.0500 | flat |
| m=400 search1500 (below n/8) | 0.0733 | 0.0733 | flat, original order |
| m=600 search1500 (above n/8) | 0.0867 | 0.0700 | 1.24× |
| m=2000 construct | 0.0500 | 0.0300 | 1.67× |
| m=2000 search1500 | 0.1400 | 0.0800 | **1.75×** |

The unconditional transpose was measured first and **regressed small m**
(m=10 search 0.0500 → 0.0700): S's m columns are the same every iteration and
stay resident, while the recomputed points' columns are cold. That is what the
`m > n/8` switch exists for, and it is why the lever ships as a switch rather
than a replacement.

**Contract.** The matrix the kernel receives from `.AsDistMatrix()` is now
exactly symmetric, since the kernel reads whichever triangle is cheaper and
that choice depends on `m`: a matrix symmetric only to rounding could
otherwise answer differently for different `k`. `SymmetryScan_cpp`
measures the largest discrepancy tile-by-tile (so a tile and its transpose are
resident together); zero passes through untouched, anything up to
`options(MaxMin.symmetryTolerance=)` has its triangles averaged with an
immediate warning, and the rest is refused. The finiteness scan runs first,
since `NA != NA` would otherwise report a missing value as asymmetry.
`MaxMean()`/`MeanDist()` pass `symmetric = FALSE`: they average the two
triangles themselves, which is documented behaviour and unaffected by this.

The k-centre solvers' own `IsSymmetric_cpp` guard went with it: it ran after
`.AsDistMatrix()`, whose guarantee is strictly the stronger one, so it could
no longer fire.

**The guard is only paid where symmetry is exploited.** It sweeps the whole
matrix — 3.0-3.5 ms at n = 2000, the single-core memory-bandwidth floor for
32 MB (12.5 GB/s measured) — which is more work than an O(nk) solve does in
total: `FarFirst(20L, d2000)` measured 6.0 ms without it and 14.0 ms with
(interleaved min-of-4, three clean `--preclean` builds). So the callers that
take `d` as written opt out with `symmetric = FALSE`: `FarFirst()`,
`PickPoint()`, `.GonzEnsemble()` (whole-column reads, and `MatrixOffDiagMax_cpp`
already scans in full rather than shortcutting to a triangle), `KCentreRadius()`,
`MaxMean()`/`MeanDist()`. It stays for DropAdd (the lever), Grasp and the
k-centre solvers, which do take a triangle — and whose kernels dominate it.
`FarFirst(20L, d2000)` returns to 6.0 ms against a 6.0 ms base; the
`"diameter"` seed to 5.0 ms against 5.5 ms.

**The two sweeps are fused.** `SymmetryScan_cpp` answers finiteness and
asymmetry from one pass, so the callers that check symmetry read the matrix
once rather than twice: intake at n = 2000 7.0 -> 5.0 ms, `Grasp(20L, d2000)`
25 -> 20 ms. It carries `AllFinite_cpp`'s bitwise exponent test into the pair
loop and scans the diagonal separately, the pair sweep stepping over it. Both
reductions are order-independent, so it threads with `.NThreads()` at no cost
to reproducibility (a full-matrix scan is 2.0 ms at one thread, 0.5 ms at
four) -- though `mc.cores` is 1 by default, so that is upside, not the
default path.

**No benchmark cell exercised the round-13 lever.** `DropAdd(20L, d500)` has
m = 20 against n/8 = 62.5, and the points cell never reaches the matrix
kernel, so CI could only ever report NSD on the transposed reads -- which it
did. A cell at m = n/2 is added: `DropAdd(250L, d500, plateau = 2000L)`,
measured 15.5 -> 12.5 ms against the pre-round-13 base.

**Verification.** The 295-case trajectory battery is **bit-identical** to the
pre-change kernel — after its asymmetric axis was replaced with a second
symmetric tie-dense shape, since asymmetric input is no longer supported
input. Run against the OLD battery the change matches on all 267 symmetric
cases and differs on 28 of the 45 asymmetric ones: the lever is exact exactly
where the contract now holds. Full suite 6304 pass / 0 fail / 2 skip.

**Leads left on the table, with their estimates now stale:**
- Lazy second-minimum record. Round 10 estimated ≤1.3× against a recompute
  branch that was ~37% of the m = n/2 profile; that branch has just become
  1.75× cheaper, so the remaining upside is smaller than the estimate and the
  exactness surface (count semantics reproducing rescanned values) is
  unchanged. Re-measure the branch share before building it.
- K-row coordinate pre-gather for the points recompute (~8% of the m = 10⁴
  profile, exact by copy). Untouched by this round.

Status: Area 2 → OPTIMISED; the matrix kernel's serial-at-limit certification
from round 10 stands for its passes (argmax fusion and threading both measured
flat there) but not for the recompute branch, which this round moved.

## Round 14 — 2026-08-18 — Area 2 (DropAdd): the last symmetry site was dead
code; removed, and the drivers moved onto the production protocol

**Trigger:** user decision. With asymmetric input refused at intake (round 13),
the matrix kernel may treat `d` as exactly symmetric everywhere, so this round
swept it for any remaining site where symmetry removes a read.

**There is exactly one, and it was API-dead.** Every other matrix read in the
kernel is a whole column that is genuinely needed: the drop and add passes
touch `d(i, x)` for all `i`, and the recompute branch's layout choice was
already taken in round 13. The one full-matrix sweep left was the
construction's max-row-sum seed (`seed0 = -1`), and `DropAdd()` had not
reached it since `0214ab2`: the wrapper fills `matrixSeed0` from
`.PickPoint(dmat, "peripheral")`, which returns >= 1 on every branch, so the
`seed0 >= 0` arm was always taken. `0214ab2` made that change deliberately —
the max-row-sum anchor was both O(n^2) *and* the worst of the seven profiled
(mean gap to the proven optimum 0.029 against the peripheral anchor's 0.011
over a 40-cell grid). Only `.DropAddTrace()` and the profiling drivers still
reached it.

**The lever was built and verified before the reach problem surfaced,** so its
cost is recorded rather than guessed. Column `j`'s entry `i < j` serves both
`rs[i]` and `rs[j]`, so the upper triangle alone suffices — `n(n+1)/2` elements
rather than `n^2` — with `rs[j]` accumulated in a register chain stored once,
complete, at the diagonal. Summation order is untouched: `rs[k]` still
receives `d(k, 0..k-1)`, then `d(k, k)`, then `d(k, k+1..n-1)`, the same
left-associated chain, exact because `d(i, j) == d(j, i)` exactly. Columns
band four at a time, since the chain is a loop-carried dependency and one
chain per column would trade the halved traffic for a serial add latency
costing more than it saves. (Round 9's banding result does not refute this:
that loop was memory-bound, where banding cannot help; here banding breaks a
latency chain.) **295/295 battery bit-identical, suite green** — then declined,
and the sweep it optimised removed instead. Never timed: the disposition was
reachability-based, so a measurement would only have priced a branch no caller
reaches.

**What shipped.**
- `src/dropadd.cpp`: the max-row-sum fallback is gone. The kernel now requires
  `seed0` in `[0, n)` and errors otherwise — a range the wrapper already
  validated, so this is a tightening for direct callers only. The points
  kernel's `-1` anti-centroid fallback is untouched: it **is** the points
  path's production default.
- `.DropAddTrace()` and `.DropAddConstruct()` seed at the peripheral anchor,
  so the trace helper and the pure-R twin walk the trajectory `DropAdd()`
  actually walks. `.DropAddConstruct()` gained a `first =` argument mirroring
  `.DropAddConstructColumn()`. The oracle test's compensating row-sum seed
  went with it.
- `dropadd-timing.R` and `dropadd-vtune10.R`: matrix cells take the production
  peripheral seed, computed once outside the timed thunks. **Every matrix
  timing from rounds 4-13 measured `seed0 = -1`** and billed each construct
  cell for an O(n^2) sweep the wrapper had stopped running, so those numbers
  are not comparable with anything measured after this round. Points cells are
  unaffected — their protocol was already production.
- `dropadd-battery.R`: the matrix cases move from `-1L` to fixed explicit
  seeds (0/1/3, keeping all three shapes); the `recP(..., -1L)` cases keep the
  default, exercising the live points fallback. Verified by capturing the
  updated script against the pre-removal kernel and comparing after: 295/295
  bit-identical, within-build invariants OK. **This commit is the new
  frozen-baseline reference** — the old `-1L` matrix cases are unreproducible
  under the new script by design.

**Re-baselined on Hamilton** (`cn059`, EPYC 7702, r/4.5.1, serial; three
whole-script reps, objectives identical throughout). Full table in
baselines.md; the structural figure is `matrix n=4e3 m=10 construct` at
**0.2 ms**. Under the old protocol that cell was dominated by the row-sum
sweep, so it was almost entirely timing a warm start `DropAdd()` had already
abandoned; what remains is the ten column passes the construction performs.
No cell here is a speedup over rounds 4-13 — the matrix cells changed
protocol and the points cells changed machine — and the log says so rather
than banking a number the code did not earn. Rep spread is at or below 1 ms
on every matrix cell, tighter than the Windows box's ±20-35%.

**Refuted by design — fusing the seed row-sums into `.AsDistMatrix`'s symmetry
scan.** The scan already reads every element, but it runs on the *pre-averaged*
matrix, so scan-time row sums would be sums of the wrong matrix whenever
intake repairs a rounding asymmetry; and it would tax Grasp and k-centre
intake for a DropAdd-only benefit that no longer exists.

Status: Area 2 → OPTIMISED, unchanged. The matrix kernel is at its symmetry
limit on every live path, and the triangle-sweep lever is closed permanently —
its target no longer exists. Round 13's leads stand: lazy second-minimum
record (re-measure the branch share first) and the K-row coordinate pre-gather
for the points recompute.

last_focus: 2
