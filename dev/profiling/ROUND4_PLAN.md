# Grasp performance with NO bit-identity constraint

**Goal:** Improve the performance of GRASP as much as you can.
The bit-identity constraint that governed rounds 1–3 is **explicitly lifted**.
GRASP is a heuristic; the user is after a good result quickly, and accepts
that selections may change.

---

## 0. Where you are

- **Branch:** `perf/grasp-unconstrained`
- **Worktree:** `C:/Users/pjjg18/GitHub/worktrees/MaxMin/grasp-perf`
- **Do NOT change the working tree of `C:/Users/pjjg18/GitHub/MaxMin`**

You are building on top of:
- **T-014**: two bit-identical hoists giving 2.03–3.79× at k=100. See
  `dev/profiling/findings.md` rows T-014/T-015 and `dev/profiling/log.md`.

---

## 1. The objective changes shape — read this before optimising anything

With bit-identity gone, wall-clock at fixed output is **no longer a valid metric**.
"5× faster" is trivially achieved by searching less. The metric becomes
**quality per unit time**.

Consequences:
- A change that only reduces work per iteration (no trajectory change) can be
  verified as before, cheaply.
- Any change that alters the search trajectory — and the main lever below does —
  must be **measured** on quality-vs-time, not asserted.
- **Build the harness in §5 before writing kernel code.** Without it you cannot
  tell an optimisation from a regression.

---

## 2. Measured cost structure

From an instrumented build, n=2000, dim=10, eliteSize=10 (the canonical shape:
GRASP runs coreset-thinned at `maxCandidates = 2000`, k ∈ {10,100}, plateau ladder
2048..1 with canonical ≈ 512):

| plateau | iters | full | phase A+C | phase B | LS min-to-set ops | PR min-to-set ops |
|---|---|---|---|---|---|---|
| 8   | 10  | 1.42 s  | 1.27 s | 10.6% | 3.41e8 (28%) | 8.72e8 (72%) |
| 64  | 146 | 4.06 s  | 1.27 s | 68.7% | 2.10e9 (76%) | 6.67e8 (24%) |
| 256 | 748 | 17.51 s | 1.28 s | 92.7% | 1.03e10 (95%) | 5.81e8 (5%)  |

Key facts, all measured, all load-bearing for the plan:
- **The two hot terms swap dominance with `plateau`.** Low plateau → phase C (path
  relinking) is ~90% of wall time. High plateau (incl. canonical) → phase B's local
  search is ~93%. Any change must be evaluated at *both* ends of the ladder.
- **`crit`/pass is exactly 2.00** — the critical positions are the two endpoints of
  the minimum-distance edge. This kills several obvious levers (see §6).
- The three O(m²) per-pass loops in the local search (`di`/`gmin`, `dstar`,
  `pair_count`) are together ~2.6% of inner ops. **Not worth touching.**
- After T-014, the local-search candidate scan is
  `crit`(2) × (n−m) × (m−1) strided reads per pass. That is the floor to beat.

---

## 3. The binding constraint is the TIE-BREAK, not bit-identity

This is the central insight of the round.

The extended-improvement rule picks, among candidates achieving the *same* new
objective `best_dstar`, the one with the fewest pairs at that distance
(`pair_count`). Ties are common — `nd == dstar` whenever `cross >= dstar`.
Because the rule must find the minimum over all ties, **every one of the (n−m)
candidates must be enumerated**. No bookkeeping avoids that while the rule stands.

So the speed ceiling is set by a *behavioural* rule, not by the arithmetic. Lifting
bit-identity is what lets you change that rule.

---

## 4. Lever A — incremental nearest-selected: O(n) per pass instead of O(n·m)

**The idea.** For each out-of-selection point `s`, maintain `min1[s]` = distance to
its nearest selected point, and `arg1[s]` = which selected vertex that is.

**Maintenance under one swap** (drop `d`, add `a`) — this is what makes it cheap:
- If `arg1[s] != d`: new `min1 = min(min1, d(s,a))`; new `arg1` is `a` if `d(s,a)`
  won, else unchanged. **Exact and O(1).** Note `d(s,a) = dptr[s + a*n]` is read
  *sequentially* in `s` — cache-friendly, unlike today's strided scan.
- If `arg1[s] == d`: the nearest was just removed, so you need the second smallest
  → rescan that `s` in O(m). **Only ~(n−m)/m points have any given vertex as their
  nearest** (≈19 of 1900 at n=2000, m=100), so the total rescan cost is O(n) per
  pass, not O(n·m).

**Screening.** For a given critical drop, a candidate can *strictly improve* only if
its min-distance to the remaining selection exceeds `best_dstar`:
- `arg1[s] != drop` → that distance is exactly `min1[s]`;
- `arg1[s] == drop` → it is the second smallest, computed on demand for that small
  subset.

Screening is then a sequential scan of a plain `double` array — **no `D()` access at
all** — instead of `(n−m) × (m−1)` strided reads.

**Expected gain:** the dominant term goes O(n·m) → O(n), up to ~m/2 on that term.

**What it costs:** you can no longer enumerate ties, so the extended-improvement
tie-break must change. **That is the trajectory change, and the whole quality risk
of the round.** See §7.1.

**Correctness note to preserve:** `min` over the remaining selection is ≥ `min1`
always, and *equals* `min1` exactly when `arg1 != drop`. Getting this backwards is
the easy bug — `min1` is a lower bound in general, so you cannot prune on it
directly without the `arg1` test.

---

## 5. Harness — build this FIRST

Two comparisons, both needed:

- **(a) Equal `plateau`** — speed at comparable search effort. Detects raw
  throughput gains.
- **(b) Equal wall-clock budget** (`maxSeconds`) — **quality at equal time.** This
  is the real frontier question and the one that decides whether the round is a win.
  A lever that wins (a) but loses (b) is a regression.

Design notes:
- Cases: canonical coreset shapes (n=2000 thinned, k ∈ {10,100}) plus real cases
  from FurthestPoint `data/cases.rda` if loadable (`../furthest-point`).
- Several seeds per cell — GRASP is randomised, so report mean **and worst** T_k,
  not just mean. A lever that raises the mean while occasionally collapsing is bad
  for a manuscript claim.
- Baseline lib = PR3 tip build; candidate lib = your build. Reuse the existing
  scaffolding: `dev/profiling/drivers/grasp-battery.R` (grid capture/compare) and
  `dev/profiling/drivers/grasp-timing.R` (clean-build timing, best-of-3).
- The old battery asserts bit-identity. **Fork it**, don't repurpose it: keep the
  identity check available for any sub-lever that *is* identity-preserving, and add
  a quality-vs-time comparison for those that aren't.

---

## 6. Already measured and rejected — do not redo (T-015)

- **The two-smallest hoist inside the local search.** Because `crit` is exactly 2
  there is only ~2× of op count to win, and the `near_two` update costs more per
  element than the branchless `if (v < best)` min it replaces. Net **+1.17× at
  k=100 but 0.64× at k=10** — a regression on a canonical cell.
- **A branchless row-buffer rewrite of the above**, pulling the strided row into an
  L1-resident buffer. **Worse still** (16.05 s vs 11.13 s shipped, plateau=256).
  The reads it targeted are already cheap: consecutive candidates share cache lines
  down each column. Do not re-attempt on cache-locality grounds.
- **Incremental maintenance of *both* minima across passes** in the T-014 sense — a
  swap can evict the argmin and repairing `min2` exactly needs the 3rd smallest and
  beyond. Note lever A in §4 sidesteps this by rescanning only the small subset
  whose `arg1` was evicted, and by not needing `min2` for the bulk of candidates.

---

## 7. Open decisions — resolve with the user, do not silently pick

1. **Tie-break replacement.** The quality risk of the round. Options: drop it (take
   the first strict improvement); keep it but only among screened winners; a cheap
   surrogate. The extended-improvement rule is a plateau-escaping device from
   Resende et al. (2010), so removing it outright may cost quality — measure with
   §5(b) before committing.
2. **`.Grasp_R` parity.** `.Grasp_R` is the pure-R spec and the parity test
   (`tests/testthat/test-grasp.R` test 2) is the best bug-catcher in the package —
   exactly when the C++ is about to get intricate. Either update it in lockstep
   (costly, keeps the oracle) or retire the test (cheaper, loses it).
   **Recommendation: keep lockstep.**
3. **FurthestPoint pin + cached figure data.** Results change deliberately, so the
   manuscript's cached GRASP figure data is invalidated. The canonical runs pin
   **MaxMin@`2a83a91`** at `../furthest-point/scripts/DEPLOY.md:16`; the user's
   planned canonical re-run must re-pin. Confirm scope before assuming.

---

## 8. Then the floor moves to construction

Once the local search is O(n)/pass, `grasp_construct` dominates: it is O(n·m) per
GRASP iteration (m greedy steps, each an O(n) argmax plus an O(n) update). At
n=2000, m=100 that is ~2e5 vs ~1.2e4 for the new local search. **So realistic
end-to-end is ~10× before construction binds** — and the 5× target likely lands
here, not beyond.

Beating construction needs a lazy max-heap greedy: `g` only ever *decreases*, so a
heap with stale-entry rejection (pop until the top's cached value matches current
`g`) is valid and can beat O(n) per step in practice. Second phase — attempt only
if lever A lands and is verified.

---

## 9. Environment gotchas — all hard-won, all will bite

- **`R` is aliased to `Invoke-History` in PowerShell.** `R CMD INSTALL ...` fails
  with a bizarre positional-parameter error. Run `R CMD *` through the **Bash**
  tool, not PowerShell.
- **Run `Rscript` through the PowerShell tool, not Bash** — MSYS bash segfaults or
  returns empty output for `Rscript` on this box. (Yes, this is the opposite of the
  previous bullet. Both are real.)
- **PowerShell escaping:** `$` inside a double-quoted `Rscript -e` string breaks.
  Put R code in a `.R` file and run the file.
- **The installed MaxMin DLL is locked by the human's running RStudio sessions**, so
  `R CMD INSTALL` into the default library fails with "cannot remove earlier
  installation" — and its restore step can leave the library broken. **Always
  install to a scratch library** via `R CMD INSTALL --no-multiarch --preclean -l
  <scratchlib> .` and set `R_LIBS_USER=<scratchlib>;<real lib>`. You need two
  scratch libs anyway (baseline vs candidate). Never kill the human's rsession
  processes.
- **Never benchmark via `devtools::load_all()`** — it compiles Rcpp at `-O0` and
  misreports speed by 8–10×. Always `R CMD INSTALL --preclean` then a fresh session.
- **`testthat::test_dir()` produces 6 spurious errors** in `test-maxentropy.R` and
  `test-print.R` — unqualified internal-function lookups, a harness artifact, not
  real failures (verified: the functions exist in the namespace). Expect
  **4893 pass / 0 fail / 6 err** as the green baseline.
- **Timing policy:** local runs give *relative* ratios only, for choosing between
  implementations. **All canonical wall-clock figures belong on Hamilton** (use the
  `/hamilton` skill) — this Windows box is not the benchmark hardware.
- **Commit identity:** per the global rules, commit as `ms609-agent` via
  `git -c user.name="ms609-agent" -c user.email="313734811+ms609-agent@users.noreply.github.com"`.
  But note `ms609-agent` is currently **blocked**, so pushes and PRs go under the
  human's own account (they instructed this explicitly for PR3).
- House style: exported fns `CamelCase`, internal helpers `.`-prefixed, never `<<-`.
  Update `NEWS.md` and bump the `.900X` dev suffix after a user-visible change.
  New R/C++ code is expected to reach 100% coverage (`coverage-gate` skill).

---

## 10. Definition of done

- Quality-vs-time harness exists and is committed under `dev/profiling/drivers/`.
- Lever A implemented, with the tie-break decision (§7.1) explicitly recorded.
- **Quality measured, not assumed**: equal-budget comparison across seeds and both
  ends of the plateau ladder, reporting mean and worst T_k.
- New `T-016`+ rows in `dev/profiling/findings.md` and a round entry in
  `dev/profiling/log.md`, in the existing house format (evidence column with real
  numbers; rejected levers recorded as NOTE rows so they are not re-attempted).
- `NEWS.md` entry stating plainly that **selections change**, and the version bumped.
- Full suite green at the 4893/0/6 baseline; `.Grasp_R` decision honoured.
- Honest reporting of the speedup achieved, even if it is short of 5×. §8 suggests
  ~10× on the local search is reachable but construction then binds.
