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
- **Remaining OPEN (out of scope this pass):** MF-03 (timeBudgetS overshoot), MF-04 (int overflow UB), FF-005/FF-006 (design/info), T7-01/T7-03/T7-04/T7-06/T7-07/T7-09/T7-10/T7-12/T7-14/T7-15/T7-16 (weak test gaps)
- R CMD check: **Status OK** (0 errors, 0 warnings, 0 notes)

### Fix pass — 0.0.0.9002  (2026-06-10)
All remaining OPEN findings from Rounds 1–7 addressed:
- FIXED (code): MF-04 (int→long long overflow UB in dropadd.cpp/dropadd_mf.cpp/grasp.cpp)
- FIXED (doc): MF-03 (timeBudgetS 256-iter overshoot documented); FF-006 (asymmetric matrix accepted, now documented)
- FIXED (behaviour): FF-005 (score=NA_real_ on all-points early return)
- FIXED (test gaps/weakness): T7-01, T7-03, T7-04, T7-06, T7-07, T7-09, T7-10, T7-12, T7-14, T7-15, T7-16
- R CMD check: Status OK (0 errors, 0 warnings, 0 notes)
- No findings remain OPEN in Rounds 1–7.

last_focus: 7
