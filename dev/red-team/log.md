# Red-team rotation log — MaxMin

Per-round notes for the `/red-team` skill. Newest entries appended at the bottom.
The `last_focus:` pointer at the very bottom tells the next invocation which area was last reviewed;
the next area is `(last_focus mod N) + 1` where N = number of areas in `focus-areas.md` (currently 7).

## Per-round entry format

```
### Round <n> — area #<N> <name>  (<YYYY-MM-DD>)
- tier: <sonnet|opus|fable>            # model the finder ran at
- finder yield: <count>                # confirmed non-trivial findings (after verification)
- verifier verdicts: <e.g. 3 candidates → 2 REAL, 1 REFUTED>
- trivial fixes: <count or list>
- seam: <"still yielding" | "ran dry">
- high-sev signal: <yes/no — detail if yes>
- next-visit decision: <"stay at <tier>" | "escalate to <tier>" | "escalate (high-sev signal)">
```

The skill reads the most recent entry for an area to set `last_tier:` and `yield:` for the escalation
decision. Escalation rules (SKILL.md §Model tiers / Normal run):
- yielded → re-visit SAME tier with a fresh agent (mine the seam)
- ran dry → escalate one tier (sonnet → opus → fable)
- high-sev signal → escalate immediately regardless

---

### Round 1 — area #1 DropAdd (single-fidelity)  (2026-06-10)
- tier: sonnet
- finder yield: 2                       # DA-01, DA-02 confirmed (after verification)
- verifier verdicts: 3 candidates → 2 REAL (DA-01, DA-02), 1 REFUTED (DA-03 — eps tie-break scaling was O(m²) not O(m³ᐟ²); docstring stands)
- trivial fixes: none
- seam: still yielding
- high-sev signal: no
- next-visit decision: stay at sonnet (seam still yielding — both findings low-severity; a fresh Sonnet pass may surface more before escalating)
- area #1 last_tier: sonnet
- area #1 yield: 2

### Rounds 2–7 — full Sonnet sweep (areas 2–7, run concurrently as an unattended batch)  (2026-06-10)
All six areas reviewed at `sonnet` (first visit). Finders reported trivial fixes rather than applying them
(parallel-safety). Verification: high/medium-confirmed claims → opus peer verifier; low/info → haiku.
Every candidate this batch verified REAL.

### Round 2 — area #2 DropAdd (multi-fidelity)
- tier: sonnet · yield: 3 (MF-01, MF-03, MF-04; all low) · verifier: haiku, 3 REAL
- trivial pending: MF-02 (want_trace coverage artifact, ← DA-01)
- seam: still yielding · high-sev: no · next: stay at sonnet
- area #2 last_tier: sonnet · yield: 3

### Round 3 — area #3 GRASP
- tier: sonnet · yield: 4 (GRASP-01 HIGH, GRASP-02 HIGH, GRASP-04 low, GRASP-05 info) · verifier: opus (01/02 — live SIGSEGV reproduced), self (04/05)
- trivial pending: GRASP-03 (dead var candMask)
- seam: still yielding · **high-sev: YES — GRASP-01/02 confirmed segfault at documented alpha=1 / alpha>1** · next: stay at sonnet (seam still yielding cheaply); the HIGH bugs are filed for fixing, not a reason to escalate the finder
- area #3 last_tier: sonnet · yield: 4

### Round 4 — area #4 FarFirst
- tier: sonnet · yield: 6 (FF-001 med-high; FF-002/003/004 low; FF-005/006 info) · verifier: opus (FF-001 reproduced), haiku (rest)
- seam: still yielding · high-sev: FF-001 med-high (silent duplicate selection) · next: stay at sonnet
- area #4 last_tier: sonnet · yield: 6

### Round 5 — area #5 Exact + maximin core
- tier: sonnet · yield: 1 (E5-01 low; E5-02 info/non-bug) · verifier: haiku, REAL
- seam: **ran dry — nothing new at this tier** (only a defensive C++-guard gap; the exact solver's binary-search / feasibility / proven-flag logic traced clean)
- high-sev: no · **next: ESCALATE to opus** (dry cheap seam on the highest-magnitude numerical core → subtle bugs, if any, are below Sonnet's cliff)
- area #5 last_tier: sonnet · yield: 1 (dry)

### Round 6 — area #6 Seed + score
- tier: sonnet · yield: 4 (F-601 low-med; F-602/604/605 low) · verifier: opus (F-601), haiku (rest)
- trivial pending: F-606 (cosmetic guard)
- seam: still yielding · high-sev: no · next: stay at sonnet
- area #6 last_tier: sonnet · yield: 4

### Round 7 — area #7 Test-suite health
- tier: sonnet · yield: 12 (T7-08 high-gap; T7-01/04/06/07/09/10/12 med; T7-03/14/15/16 low) · verifier: haiku, 17/17 REAL (incl. trivial)
- trivial pending: T7-02/T7-05 (redundant assertions), T7-11/T7-13/T7-17 (minor)
- seam: still yielding · high-sev: T7-08 (= FF-001 coverage gap) · next: stay at sonnet
- area #7 last_tier: sonnet · yield: 12

---

### Fix pass — 0.0.0.9001  (2026-06-10)
All P0/P1/P2 findings from Rounds 1–7 that were in scope addressed in a single batch commit:
- **FIXED (P0 — crash/corruption):** GRASP-01, GRASP-02 (segfault at alpha=1/alpha>1); FF-001 (NA dist matrix silent dup); MF-01 (NaN DropAdd_points_cpp guard)
- **FIXED (P1 — API hardening):** FF-002, FF-003, FF-004 (method/m validation); F-601 (misspelled anchor silently dropped); F-602 (winning_strategy all-NA); F-604, F-605 (MinDist NA/dup idx)
- **FIXED (P2 — defensive guards):** E5-01 (MaximinFrom_cpp/MaximinFromPoints_cpp n range); GRASP-04 (m=1 loop); GRASP-05 (alpha=1 parity test gap); T7-08 (NA matrix test gap)
- **FIXED (cosmetic/trivial):** GRASP-03 (dead candMask), F-606 (n==0L alignment)
- **Remaining OPEN (out of scope this pass):** MF-03 (maxSeconds overshoot), MF-04 (int overflow UB), FF-005/FF-006 (design/info), T7-01/T7-03/T7-04/T7-06/T7-07/T7-09/T7-10/T7-12/T7-14/T7-15/T7-16 (weak test gaps)
- R CMD check: **Status OK** (0 errors, 0 warnings, 0 notes)

### Fix pass — 0.0.0.9002  (2026-06-10)
All remaining OPEN findings from Rounds 1–7 addressed:
- FIXED (code): MF-04 (int→long long overflow UB in dropadd.cpp/dropadd_mf.cpp/grasp.cpp)
- FIXED (doc): MF-03 (maxSeconds 256-iter overshoot documented); FF-006 (asymmetric matrix accepted, now documented)
- FIXED (behaviour): FF-005 (score=NA_real_ on all-points early return)
- FIXED (test gaps/weakness): T7-01, T7-03, T7-04, T7-06, T7-07, T7-09, T7-10, T7-12, T7-14, T7-15, T7-16
- R CMD check: Status OK (0 errors, 0 warnings, 0 notes)
- No findings remain OPEN in Rounds 1–7.

last_focus: 7

---

### Round 8 — area #7 k-centre (CDSh + exact)  (2026-06-11)
- tier: sonnet (finder) → opus (KC-001 high-sev verify), haiku (KC-002..007 batch verify)
- finder yield: 7 confirmed (all REAL after verification)
- New code (added this session): `src/kcentre_cdsh.cpp`, `R/kcentre.R`, `tests/testthat/test-kcentre.R`.

**KC-001 (HIGH, opus-confirmed) — CDSh binary search on a non-monotone predicate.**
`KCentreCDSh_cpp`'s `achieved <= cand[mid]` feasibility test is NOT monotone in r
(the kernel's "feasible region is a suffix" comment was false), so the O(log n)
sampled candidates can skip the radius that yields the best construction. Opus
reproduced (bit-identical R port of the kernel): binary search loses to a same-seed
full scan in 196/960 (20.4%) random instances (n∈{10,15,50}), up to 48.2% over
optimum; and **loses to the Gonzalez 2-approximation in 11/960 cases**, violating
the package's own documented "at least as tight as Gonzalez" contract. Window-scan
mitigation does NOT work (misses are non-local); `nstart>1` helps but doesn't fix;
only a full scan reliably fixes it (O(n⁴), affordable only at small n).

**Other findings (all REAL):** KC-002 (med) asymmetric `d` → silent wrong answer
(upper-tri-only candidates + symmetry-exploiting kernel/IP, no guard); KC-003 (low)
`KCentreRadius(d, 0L)` returns `-Inf` silently on the matrix path (no upper-bound
idx check); KC-004 (low) ExactKCentre null-warm-start fallback reports `cand[hi]`
not the witness radius; KC-005/006/007 (coverage) no tests for the timeout/
`proven=FALSE` path, asymmetric input, or the dedup (`length<k`) case.

**Fixes applied this round (all confirmed findings):**
- KC-001: (a) **Gonzalez floor** — `KCentre` now returns the better of CDSh and a
  `FarFirst` (peripheral 2-approx) pass, *guaranteeing* the documented ≥ Gonzalez
  contract everywhere at O(nk) cost; (b) **exhaustive candidate scan when n ≤ 150**
  (full CDS over all radii, affordable at small n where the losses are worst),
  binary search above that; (c) corrected the false "suffix"/monotone comments.
- KC-002: `isSymmetric(d)` guard in `KCentre`/`ExactKCentre` (error on asymmetric;
  k-centre is a metric problem). `KCentreRadius` left asymmetric-tolerant (correct).
- KC-003: index-range check added to `KCentreRadius`.
- KC-004: `ExactKCentre` reports `radius = KCentreRadius(d, indices)` (always exact;
  fixes the fallback and inconclusive-break inconsistency).
- KC-005/006/007: tests added.

Verified: full `test_local` green (kcentre 118 assertions incl. the new regressions;
all other areas unaffected). 0/120 previously-losing instances now worse than
Gonzalez (was 11/960); CDSh/optimum mean 1.022, max 1.349 (near-optimal via the
small-n exhaustive scan). Perf at n=2004 k=20: 307 ms (T-010 kernel) → ~440 ms
end-to-end — the symmetry guard (in-place `IsSymmetric_cpp`, ~70 ms) plus the
Gonzalez floor (O(nk)) are the correctness cost; still ~3× faster than the
pre-T-010 1275 ms, and `dist` input skips the symmetry scan. The exhaustive scan
does not trigger at n=2004 (binary search + floor).

- seam status: still yielding (7 found) — next visit of this area stays at sonnet
  (fresh angle) before escalating.
- escalation decision: KC-001 high-sev was confirmed by opus and fixed; no further
  immediate escalation needed. Next rotation advances to area #8 (test-suite health).

last_focus: 7

---

## Round 9 — area #9 MaxMean (max-mean dispersion / RLTS) — opus

date: 2026-06-16
tier: opus (new area; implementer had already self-found+fixed 3 bugs: RNG not
seeded → set.seed ignored, aspiration ignoring non-tabu moves, tabu carryover
across restarts). Ran opus directly as peer to that self-review.

Finder yield: 7 candidates (no high-severity / wrong-answer bug). The objective
math was verified exact: incremental P-array vs scratch over millions of flips,
score==MeanDist on symmetric input, brute-force optima matched (RL on/off),
monotonicity, |S|>=2 guarantee, Sa-compaction fuzz (200k traces), index-overflow
audit (all *n casts to size_t), RNGScope consumes the stream once.

Verification: behavioural claims (MM-01 empty result, MM-03 score≠MeanDist)
confirmed by direct empirical probe rather than a verifier agent (a probe is
stronger evidence). Code+paper-fidelity claims (MM-02 overflow, MM-04/05/06/07 RL
init) confirmed by re-reading the kernel against the paper's Eqs. 4-10.

All FIXED 9003 except MM-07 (WONTFIX — paper Eq.5/Eq.7 conflict, code follows Eq.5,
documented in a comment):
- MM-01 (med): do-while restart loop (budget checked at END, matching paper Alg.1)
  → ≥1 restart always completes; tiny-budget result is now a classed |S|>=2 set.
- MM-02 (low): `iters` kept numeric (double), not as.integer → no NA past 2^31.
- MM-03 (med): root cause deeper than finder thought — the incremental sum_pairs
  adds the row-sum p_u, so on ASYMMETRIC d the objective drifts order-dependently
  (finder's fuzz was symmetric, missed it). Fix: symmetrize d once at the R
  boundary ((d+t(d))*0.5, bit-identical no-op for symmetric input); MeanDist
  symmetrizes to match. Now well-defined and consistent.
- MM-04/05/06 (low, RL-init fidelity): true max in Eq.4 (terminal=0); R init 0.0
  not 1.0; reward magnitude increments AFTER use (first improvement → 1, not 2).

Coverage: maxmean.cpp 149/149, maxmean.R 41/41, MeanDist 7/7, MaxMean print
methods 18/18 — all 100%. Full test_local green (63 maxmean assertions incl. the
MM-01/02/03 regressions; all other areas unaffected).

- seam status: ran dry on wrong-answer bugs (objective core verified exact); the
  7 findings were a contract bug + an overflow + an asymmetric-input wart + 4
  RL-init fidelity nits, all addressed. A fresh-angle re-visit could probe the RL
  convergence-vs-paper question (A/B on a published MDPI instance) but that is
  efficiency, not correctness. Next visit escalates to fable only if a
  correctness concern resurfaces.
- escalation decision: no high-sev signal. Rotation advances to area #1
  (DropAdd single-fidelity) next.

last_focus: 9

### Round 9 addendum — benchmark validation (2026-06-17)

A/B vs published best-known (Lai & Hao 2016 / Nijimbere et al. 2020 Table 1) on the
MDPI Type-I n=500 instances (downloaded from grafo.etsii.urjc.es/optsicom/edp/),
useRL=TRUE, 30 s/instance, set.seed(1) — see dev/ab-mdpi.R:

    MDPI1_500  best 81.277044  MaxMean 81.277044  MATCH
    MDPI2_500  best 78.610216  MaxMean 78.610216  MATCH
    MDPI3_500  best 76.300787  MaxMean 76.300787  MATCH  (all gaps < 5e-7)

Confirms the RL-initialisation fidelity fixes (MM-04/05/06) and the implementation
as a whole reach the paper's optima. Instances not committed (~64 MB); ab-mdpi.R
documents how to re-download.
