# Grasp round 5 — to the ceiling, then across the cores

**Mandate (user):** the project compares heuristics across a time/quality
trade-off, so each implementation must be absolutely optimised. Keep
optimising until the ceiling. Parallelisation is authorised, with core count
controlled per the `../TreeDist` house convention.

**Starting point:** round-4 tip (`8ae025c` on `perf/grasp-unconstrained`),
0.0.0.9007. Phase split at n = 2000, k = 100 (instrumented, seed-42 instance):

| plateau | construct | local search | path relink | total |
|---|---|---|---|---|
| 8 | 0.015 (9%) | 0.084 (52%) | 0.064 (39%) | 0.163 s |
| 256 | 0.351 (20%) | 1.335 (77%) | 0.044 (2.5%) | 1.736 s |

---

## Order of work — the RNG change is a one-way door

T-018/T-019/T-020 are **exact** (bit-identical to 0.0.0.9007) and land first,
each as its own battery-gated commit. T-021 changes the trajectory once and
**invalidates the battery baseline**; it lands strictly last, and nothing
exact follows it in this round without first capturing the new baseline.

### T-018 — one incremental structure for the whole LS header (exact)

A NearTwo-per-selected-member summary over the selection —
`(min1, min2, arg1, arg2, cnt)` where `cnt` counts partners at exactly
`min1` — serves everything the per-pass m²/2 pair sweep currently provides:

- `di[v]` = min1[v]; `dstar` = min over members of min1; witness = argmin
  member + its arg1;
- `pair_count` = Σ over members at dstar of cnt, halved (each dstar-pair is
  counted by both endpoints — exact integer);
- post-drop rescores for witness drops: `near_excl(nt[v], w)` per member —
  O(m) per query, replacing the lazy m²/2 exclusion pass.

Maintenance per accepted swap (drop D, add A), for each surviving member v:
- `arg1[v]` or `arg2[v]` == D → full rescan + recount (rare; witness
  invariant, as T-016/T-017);
- else read `d(v,D)`: **if it equals min1[v], decrement cnt** (a tied
  non-witness partner left — the easy-to-miss branch); read `d(v,A)`:
  `< min1` → shift + cnt=1; `== min1` → ++cnt; `< min2` → update min2/arg2.
- A itself gets a fresh scan.

Per-pass header work falls from ~2×m²/2 reads to O(m); per-swap cost is O(m)
plus rare rescans. The `count_pairs_within` memo-miss term (~2 × m²/2 per
pass) survives and is expected to be the next floor — if it binds, §7.1 of
round 4 (tie-break surrogate) revives as a *measured* question for the user,
via the frontier harness; it does not go into a commit.

Gates: battery 1458/1458; suite-level branch-counter check (the round-4
method) with special attention to the tie-decrement branch; interleaved
min-of-N timing from the first keep/reject call.

### T-019 — construction pass merge (exact)

Fold the argmax/gmin scan into the g-update pass (first-index tie semantics
preserved by ascending iteration). Construction goes from ~3 to ~2 O(n)
passes per step. Battery + interleaved timing.

### T-020 — LS init g-handoff (measure first)

Post-T-018, measure the share of the LS call spent building min1/arg1 (and
the new member summaries). Only if it binds, choose BY MEASUREMENT between:
value-only min1 screen + drop-column equality probe (construction untouched),
or argmin tracking inside construction's update loop. May be dropped.

### T-021 — nCores-invariant parallel GRASP (trajectory change; lands last)

**Determinism contract: the result depends on the seed only — identical at
every core count.** Core count affects wall-clock alone. This is the
strongest property the manuscript comparison can have, and it is stronger
than "reproducible given seed + nCores".

Design:
- **RNG.** Constructions stop drawing `R_unif_index` on demand. Before each
  batch, the MAIN thread pre-draws m uniforms per iteration from R's stream
  (`unif_rand`); workers consume them as `index = floor(u * size)`.
  Deliberate deviation: floor-index instead of R's rejection-sampled unbiased
  index — bias is negligible (size ≤ ~2000 against 53 bits) and it makes the
  pure-R mirror trivial (`runif(m)` + `floor`). Fixed consumption of m draws
  per construction (today's consumption is data-dependent) — documented.
- **Phase B** runs batch-synchronously with a **fixed constant B = 32**
  (part of the algorithm definition, NOT a tuning knob — R RNG consumption
  depends on it). Construct + local search are elite-set-independent, so a
  batch is embarrassingly parallel; results merge on the main thread **in
  iteration order** through `grasp_try_insert`, with plateau/maxIter/budget
  evaluated at merge time. Overshoot waste ≤ B−1 iterations per run.
  Merged-prefix results are therefore invariant to thread scheduling AND to
  B — but B stays fixed anyway for RNG-consumption stability.
- **Phase A** is the same batch treatment (eliteSize constructions).
- **Phase C** is deterministic given the elite set: the 45 pair-relinks
  parallelise directly with a reduce in pair order (first-best-wins ties
  preserved).
- **Threading:** `getOption("mc.cores", 1L)` → `n_threads` argument
  (Grasp_cpp signature changes → `compileAttributes()` + `ccache -C` before
  rebuilding); `Makevars` + `Makevars.win` with `$(SHLIB_OPENMP_CXXFLAGS)` in
  PKG_CXXFLAGS and PKG_LIBS; `#ifdef _OPENMP` guards with a serial fallback;
  per-thread scratch after TreeDist's MatchScratch pattern; no R API off the
  main thread (`cb`, `checkUserInterrupt`, RNG all between batches).
- **`.Grasp_R`** updated in lockstep (same batching, `runif`-based draws) so
  the parity test remains the structural oracle.

Gates: (i) invariance battery — one seed, nCores ∈ {1, 2, max}, bit-identical
results; (ii) parity vs updated `.Grasp_R`; (iii) frontier equal-budget
non-inferiority vs 0.0.0.9007 across seeds (the trajectory change must not
cost quality); (iv) cores-vs-speedup curve (extend grasp-frontier.R with an
nCores dimension — do not fork it); (v) suite green with tests capped at
2 cores (`_R_CHECK_LIMIT_CORES_`). Then capture the NEW battery baseline.

**Cost, stated plainly:** T-021 un-does round 4's headline. Selections change
once, so the manuscript's cached GRASP figure data is invalidated and the
FurthestPoint re-pin (round-4 §7.3) becomes owed after all. The mandate
authorises this; the payoff is nCores-invariant determinism. Goes in NEWS
(0.0.0.9008) and the round record as an explicit trade.

---

## Ceiling statement

After T-018/T-019 the serial residue is AT-LIMIT-shaped: an O(n) screen the
tie-break provably requires, Ω(n·m) construction bound by its own sequential
update pass, and the tie-break's own memo cost. From there the remaining axis
is core count. Local timings choose between implementations only; real
scaling curves and canonical wall-clock belong on Hamilton, plus a Linux/gcc
sanity check before any mission-wide claim (per the /profile skill).
