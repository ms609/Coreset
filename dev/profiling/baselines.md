# Baselines — MaxMin  (for `/profile regress`)

Median wall time, `bench::mark`, R-devel, `-O2`. Refresh each round.

## Area 1 — FarFirst (default ensemble, 3 starts, N=6000) — AFTER T-001/T-002

| dim | n | n/N | points ms | matrix ms |
|----:|--:|----:|----------:|----------:|
| 2 | 60 | 0.01 | 3.5 | 2.7 |
| 2 | 300 | 0.05 | 14.0 | 10.1 |
| 2 | 1200 | 0.20 | 52.7 | 36.9 |
| 2 | 3000 | 0.50 | 128.5 | 90.0 |
| 10 | 60 | 0.01 | 7.8 | 2.7 |
| 10 | 300 | 0.05 | 34.2 | 10.2 |
| 10 | 1200 | 0.20 | 129.9 | 37.1 |
| 10 | 3000 | 0.50 | 322.1 | 92.0 |

(Pre-optimisation figures and speedups are in log.md, Round 1.)

## Area 2 — DropAdd — AFTER T-005a (matrix seed reorder)

> **Every matrix row in this section and the two below measures the retired
> `seed0 = -1` protocol.** The cells passed the kernel's max-row-sum warm
> start, which `DropAdd()` stopped using at `0214ab2` and round 14 removed, so
> each `construct` row is billed for an O(n^2) sweep no caller ran. They are
> kept as the record of what was measured, and are **not comparable** with
> anything timed after round 14. Points rows are unaffected — that path's
> `-1` anti-centroid seed is its production default. Regress against the
> round-14 section.


| case | metric | ms |
|------|--------|---:|
| matrix construction, n=4000, m=10 (maxIter=0) | median | 10.3 |
| matrix construction, n=6000, m=10 (maxIter=0) | median | 22.4 |
| matrix search, n=4000, per-iter | — | ~0.06 |
| points search, n=20000 dim=2, per-iter | — | ~0.23 |

(Pre-fix construction: 116 ms / 293 ms → 11.3× / 13.1×.)

### AFTER round 10 (blocked fills + fused-pass mc.cores, 2026-08-13)

All trajectories bit-identical (295-case battery incl. drop/add sequences;
nCores-invariant). Kernel-direct cells, `drivers/dropadd-timing.R`,
interleaved minima; regress against THESE rows. Threads engage on the
points kernel only at n ≥ 16384 (the matrix kernel is serial by
measurement — see log.md round 10).

| case | 1 thread ms | 8 threads ms |
|------|---:|---:|
| matrix n=4e3 m=10 construct (maxIter=0) | 8 | — |
| matrix n=4e3 m=2000 construct (maxIter=0) | 45 | — |
| matrix n=4e3 m=10 search1500 | 50 | — |
| matrix n=4e3 m=2000 search1500 | 130 | — |
| points n=2e4 d2 m=10 search1000 | 217 | ~180 |
| points n=2e4 d10 m=10 search600 | 277 | ~210 |
| points n=2e4 d10 m=1e4 construct | 1460 | 610 |
| points n=2e4 d10 m=1e4 search600 total | 1670 | — |

(Pre-round-10 points cells: 233 / 327 / 1860 / 2110 ms. The 8-thread
figures are from the same-code experiment build and the public-API run;
mc.cores end-to-end at n=2e4 m=2000 plateau-capped: 320 → 180 ms, 1.75×.)

## Area 2 — DropAdd matrix kernel — AFTER round 13 (symmetric-read recompute)

Kernel-direct, n = 4000, interleaved min-of-3. Objectives identical per cell.

| cell | before s | after s |
|------|---------:|--------:|
| m=10 construct | 0.0100 | 0.0100 |
| m=10 search1500 | 0.0500 | 0.0500 |
| m=400 search1500 | 0.0733 | 0.0733 |
| m=600 search1500 | 0.0867 | 0.0700 |
| m=2000 construct | 0.0500 | 0.0300 |
| m=2000 search1500 | 0.1400 | 0.0800 |

Recompute reads switch to column-major at m > n/8. 295-case trajectory battery
bit-identical (asymmetric axis retired with the symmetry contract).

## Area 2 — DropAdd — AFTER round 14 (production seed protocol) — REGRESS HERE

Hamilton8, r/4.5.1, gcc/12.2, `OMP_NUM_THREADS=1`, `mc.cores` unset.
`drivers/dropadd-timing.R`, three whole-script reps per job, each cell an
internal min-of-3. Objectives were identical on every cell of every rep in
both jobs — they are the trajectory identity probe, so a changed objective
means a changed trajectory, not a faster one.

Two jobs are shown because they landed on different nodes of `-p shared`, and
the difference is the point (see the variance note below).

| cell | cn059 s | cn058 s | node Δ |
|------|--------:|--------:|-------:|
| matrix n=4e3 m=10 construct | 0.0002 | 0.0002 | 0% |
| matrix n=4e3 m=2000 construct | 0.0245 | 0.0247 | +1% |
| matrix n=4e3 m=10 search1500 | 0.0672 | 0.0670 | 0% |
| matrix n=4e3 m=2000 search1500 | 0.0893 | 0.0927 | +4% |
| matrix n=4e3 m=400 search1500 | 0.0620 | 0.0923 | **+49%** |
| matrix n=4e3 m=600 search1500 | 0.0820 | 0.0967 | **+18%** |
| points n=2e4 d2 m=10 search1000 | 0.3680 | 0.3680 | 0% |
| points n=2e4 d10 m=10 search600 | 0.3987 | 0.4060 | +2% |
| points n=2e4 d10 m=1e4 construct | 2.0170 | 2.0230 | 0% |
| points n=2e4 d10 m=1e4 search600 | 2.2990 | 2.2930 | 0% |

### A kernel row is not the cost of a call

Measured on cn058, same job:

| | s | share of a `DropAdd(20, d4)` call |
|---|---:|---:|
| `.AsDistMatrix(d4)` intake scan (O(n²)) | 0.0325 | **78%** |
| `.PickPoint(d4, "peripheral")` seed (O(n)) | 0.0001 | <1% |
| END-TO-END `DropAdd(20, d4, plateau=200)` | 0.0417 | — |
| END-TO-END `DropAdd(600, d4, plateau=200)` | 0.0483 | — |

The intake scan is **160× the `m=10 construct` kernel cell** and about
four-fifths of a whole small-`k` call. Round 13 kept that scan for DropAdd
deliberately — the recompute branch's triangle choice needs exact symmetry —
so it is real, required, per-call work that no kernel row contains.

`farfirst-timing.R` times `FarFirst()` through the **public API**, so its rows
include an intake of the same kind. Only the `END-TO-END` rows here compare
with it. Reading a kernel row against a FarFirst row understates DropAdd by
roughly the intake, which is most of the call at small `k`.

### Variance: interleave within one job, and distrust the switch pair

Within a job, rep spread reaches 4.0 ms (6.5% on `m=400 search1500`). **Across
jobs it is far worse**: `m=400` moved 49% and `m=600` 18% between cn059 and
cn058, on identical code walking identical trajectories. Those two cells sit
either side of the `m = n/8` recompute switch and are the most
memory-bandwidth-sensitive in the set, so on `-p shared` — where another job
can be competing for bandwidth on the same node — their absolute values are
not trustworthy to better than ~50%.

Consequence for future rounds: **do not regress a matrix cell against an
absolute number from a previous job.** Run both arms interleaved inside one
job on one node, as `coreset-kc/kcab.sh` does. The stable cells (everything
except the switch pair) hold to ~4% across nodes and can be read absolutely.

**No row here is a speedup over the sections above.** The matrix cells changed
protocol (peripheral seed, not `seed0 = -1`) so they walk different
trajectories, and the points cells changed machine. Ratios across that
boundary are meaningless in both directions.

The one figure worth reading structurally: `m=10 construct` at **0.2 ms**. The
same cell under the old protocol was dominated by the O(n²) max-row-sum
sweep, which no caller ran — the cell was almost entirely measuring a warm
start `DropAdd()` had already abandoned. What remains is the ten column
passes the construction actually performs.

## Area 3 — Grasp — AFTER T-007 (base_z min-edge witness hoist)

| case | metric | ms |
|------|--------|---:|
| n=200, m=50, eliteSize=5, plateau=15 | median | 16.4 |
| n=200, m=100, eliteSize=5, plateau=15 | median | 39.5 |

(T-006→T-007: 22.2→16.4 ms / 85.0→39.5 ms → 1.34× / 2.07×, 416/416 bit-identical.
 Cumulative vs pre-T-006: 315→16.4 ms / 1329→39.5 ms → 19.2× / 33.6×.)

### AFTER round 4 (T-016/T-016b/T-017, 2026-08-13) — canonical coreset shapes

`drivers/grasp-timing.R` best-of-3 minima across interleaved invocations
(this box spikes ±20–35% on sub-second cells — regress-compare against
interleaved minima, never one run; see log.md round 4).

| case | metric | ms |
|------|--------|---:|
| n=2000 dim=10 k=100 plateau=8 | min | 120 |
| n=2000 dim=10 k=100 plateau=64 | min | 520 |
| n=2000 dim=10 k=100 plateau=256 | min | 2390 |
| n=2000 dim=10 k=10 plateau=64 | min | 50 |
| n=500 dim=10 k=50 plateau=64 | min | 40 |

(Round 4 vs PR3 tip: 470→120 / 2210→520 / 10490→2390 ms → 3.9× / 4.3× / 4.4×,
 1458/1458 bit-identical, equal-plateau objectives identical at n = 2000.
 Quality-vs-time harness: `drivers/grasp-frontier.R`.)

### AFTER round 5 (T-018/T-019/T-021, 2026-08-13) — single-threaded

**Trajectory changed once at T-021** (batched floor-index draws): `iters`
per plateau differ from all earlier rows, so cross-version wall-clock is not
per-iteration comparable — per-iteration cost is unchanged (~2.1 ms at
plateau 256). Regress against THESE rows, single-threaded, interleaved
minima. Determinism is nCores-invariant (`drivers/grasp-invariance.R`);
scaling curve in `drivers/grasp-scaling.R` (plateau 256:
0.83/0.43/0.23/0.14 s at 1/2/4/8 threads, 16-core box).

| case | metric | ms |
|------|--------|---:|
| n=2000 dim=10 k=100 plateau=8 (iters 8) | min | 120 |
| n=2000 dim=10 k=100 plateau=64 (iters 64) | min | 220 |
| n=2000 dim=10 k=100 plateau=256 (iters 372) | min | 820 |
| n=2000 dim=10 k=10 plateau=64 (iters 64) | min | 20 |
| n=500 dim=10 k=50 plateau=64 (iters 163) | min | 50 |

### AFTER round 6 (tie-arm derivation + g-handoff, 2026-08-13) — single-threaded

Bit-identical to the round-5 rows (same iters, objectives, pr_calls;
battery 1458/1458). Regress against THESE rows, single-threaded,
interleaved minima. The k=10 row is per 40 calls (single calls sit under
this box's timer granularity).

| case | metric | ms |
|------|--------|---:|
| n=2000 dim=10 k=100 plateau=8 (iters 8) | min | 110 |
| n=2000 dim=10 k=100 plateau=64 (iters 64) | min | 150 |
| n=2000 dim=10 k=100 plateau=256 (iters 372) | min | 610 |
| n=2000 dim=10 k=10 plateau=64, 40 calls | min | 1020 |
| n=500 dim=10 k=50 plateau=64 (iters 163) | min | 40 |

Re-verified 2026-08-13 after the rebase onto post-#2 main (round 8, no code
change shipped): interleaved minima 100/150/590/—/40 ms — within noise of
the rows above. The 40-call k=10 row is seed-protocol-sensitive (fixed-seed
calls; a seeds-1:40 mix measures ~1.6 s on the same build) — regress it
only with the original fixed-seed protocol.

## Area 1 — FarFirst single pass — AFTER T-004 (column reorder, AT-LIMIT)

| case | metric | ms |
|------|--------|---:|
| points, dim=2, N=6000, n=3000 | median | 40.4 |
| points, dim=10, N=6000, n=3000 | median | 108.8 |

### AFTER round 7 (argmax fold + mc.cores kernels, 2026-08-13)

Bit-identical selections (925-case battery, both paths, every strategy).
Interleaved min-of-3, `drivers/farfirst-timing.R`; regress against THESE
rows. Threads engage on the pass only at N ≥ 32768.

| case | 1 thread ms | 8 threads ms |
|------|---:|---:|
| matrix, dim=10, N=6000, k=3000 | 233 | 235 (below threshold) |
| points, dim=2, N=6000, k=3000 | 23 | 23 (below threshold) |
| points, dim=10, N=6000, k=3000 | 87 | 90 (below threshold) |
| points, dim=2, N=100000, k=1000 | 140 | 57 |
| points, dim=10, N=100000, k=1000 | 560 | 125 |
| anchors ensemble (diameter+anti_medoid+rownorm), N=6000, k=300 | 740 | 140 |

### AFTER round 9 (validation scan + block sweeps + anchor scans, 2026-08-13)

Bit-identical selections and scores everywhere (battery 925/925; cross-path
identity at dims 1–11; nCores-invariant). Interleaved min-of-3,
`drivers/farfirst-timing.R` (which now includes the matrix-anchors cell);
regress against THESE rows. The matrix cells' validation scan
(`AllFinite_cpp`) follows mc.cores; the matrix greedy pass itself is serial
at every RAM-feasible N (round 9 lever F).

| case | 1 thread ms | 8 threads ms |
|------|---:|---:|
| matrix, dim=10, N=6000, k=3000 | 42 | 28 |
| points, dim=2, N=6000, k=3000 | 18 | 20 |
| points, dim=10, N=6000, k=3000 | 67 | 67 |
| points, dim=2, N=100000, k=1000 | 107 | 47 |
| points, dim=10, N=100000, k=1000 | 395 | 85 |
| anchors ensemble (diameter+anti_medoid+rownorm), N=6000, k=300 | 280 | 60 |
| matrix anchors ensemble (same anchors), N=6000, k=300 | 200 | 90 |

### AFTER round 10 (restart-parallel ensembles, 2026-08-14)

Bit-identical selections and scores everywhere (battery 925/925 at mc.cores
1/2/8; nCores-invariant). The rows above are unchanged — this round moved
only the multi-start cells, which the timing driver gains here (`ens-default*`,
seeded per call: the random-furthest draw must not vary between A/B arms).
Interleaved min-of-5, rep counts sized so each timed block clears this box's
~10-16 ms clock granularity; regress against THESE rows.

| case | 1 thread ms | 4 threads ms | 8 threads ms |
|------|---:|---:|---:|
| ens-default matrix, N=6e3, k=3000, nSeeds=3 | 113 | 42 | 40 |
| ens-default points, dim=10, N=6e3, k=3000, nSeeds=3 | 205 | 72 | 72 |
| ens-default matrix, N=6e3, k=3000, nSeeds=8 | 206 | 66 | 50 |
| ens-default points, dim=10, N=6e3, k=3000, nSeeds=8 | 540 | 153 | 123 |

(Round 9 → 10 at 8 threads: 80→40 / 207→72 / 174→50 / 547→123 ms —
**2.00× / 2.88× / 3.48× / 4.43×**. Serial within noise: 117→113 / 205→205 /
206→206 / 553→540. The anchors ensembles are unmoved by design — their
O(N²) seeds dominate and already parallelise: matrix anchors 106→106 ms.)

## Area 4 — ExactMaxMin — AFTER T-008 (sparse-A + Grasp warm-start gallop)

Wall time per solve (single call, highs threads=1, `set.seed(1)`); OLD = dense-A
full-bisection solver. Whole-grid totals over 14 ground-truth cases × k∈{2,4,6,10}.

| case (k=10) | n | OLD s | NEW s | speedup |
|------|--:|------:|------:|--------:|
| tc4_hierarchical | 240 | 1.03 | 0.10 | 10.3× |
| tc22_penguins | 342 | 7.81 | 0.50 | 15.6× |
| tc11_ionosphere | 351 | 16.78 | 0.97 | 17.3× |

Grid total (56 solves): OLD 222.6 s → NEW 11.3 s (**19.6×**); per-case median 16.3×.
k=2 solves need **zero** IP solves (heuristic attains the diameter) → ~100–1300×.
Correctness: 56/56 proven optima bit-identical to OLD; brute-force oracle + full
suite (598) green. Worst case tc20_zoo k=4 (n=101): 0.59→0.67 s (0.88×, warm-start
overhead on a tiny instance; irrelevant at the n≥342 job sizes).

### AFTER round 11 (combinatorial probe reduction + capped clique IPs)

Base = perf/dropadd fb93c2d; R-devel 4.7.0; highs 1.14.0.2 (scratch lib);
interleaved min-of-5 via `drivers/exact-grid.R` (seed 1000+cell). Scores and
proven flags identical 56/56; witnesses free to differ (user contract).

| measure | fb93c2d | round 11 | speedup |
|------|------:|------:|--------:|
| grid total (56 cells) | 10.40 s | 5.42 s | 1.92× |
| tc11_ionosphere k=10 | 920 ms | 170 ms | 5.4× |
| tc22_penguins k=4 | 680 ms | 340 ms | 2.0× |
| tc10_glass k=10 | 180 ms | 10 ms | 18× |
| tc8_highdim_gaussians k=6 | 670 ms | 500 ms | 1.34× |
| tc1_uniform k=10 (known trade) | 140 ms | 270 ms | 0.52× |

Stretch (single pass, 600 s budget, shipping build with setup lever;
interleaved min-of-3 for the spam/satellite comparisons): breastcancer k10
4.0 s, pima k10 3.6 s, vehicle k10 20.6 s, vowel k10 172 s; spam n=4601
k4/10 = 1.8/2.0 s (fb93c2d: 201 s and 8.9 GB peak at k=4); satellite
n=6435 k=4 = 3.9 s (fb93c2d model ≈16 GB: not attempted). Spam and
satellite decide every probe combinatorially — no IP is ever posed.

## Area 4 — ExactMaxMin — AFTER round 12 (clique search replaces the IP)

Grid = 14 ground-truth cases × k∈{2,4,6,10}; interleaved min-of-3. Stretch
cells are single runs. Base = perf/exact tip fa482b3 (round 11).

| cell | base s | new s |
|------|-------:|------:|
| grid, 56 cells | 6.5 | 0.9 |
| breastcancer k10 (n=683) | 4.0 | 0.1 |
| pima k10 (n=768) | 3.6 | 0.1 |
| vehicle k10 (n=846) | 20.4 | 0.1 |
| vowel k10 (n=990) | 171.4 | 0.2 |
| spam k4 (n=4601) | 2.3 | 0.3 |
| spam k10 | 2.4 | 0.5 |
| satellite k4 (n=6435) | 4.5 | 1.0 |
| satellite k10 | — | 13.1 |

56/56 scores bit-identical to base, `proven` identical, witnesses valid. No
`highs` and no `Matrix` on the solve path. Root-branch parallelism measured
(3.36× at 8T on infeasibility proofs, flat on feasible probes, ~1.2×
end-to-end) and declined — see log.md round 12.

## Area 5 — KCentre (CDSh) — AFTER T-010 (cache reorder + C++ candidates)

Clustered Gaussian n=2004, dim=10, k=20; per-call median (6 reps), `-O2`.

| case | OLD ms | NEW ms | speedup |
|------|-------:|-------:|--------:|
| KCentre, n=2004, dim=10, k=20 | 1275 | 307 | 4.16× |

Centres + radius bit-identical at k∈{5,20,50}. `test_local(filter=kcentre)` 70/70.
ExactKCentre baseline pending (T-011).

**Update (red-team Round 8, KC-001 fix):** `KCentre` end-to-end is now ~440 ms at
n=2004 k=20 (matrix input), not 307 ms — the added correctness guards (in-place
`IsSymmetric_cpp` ~70 ms; Gonzalez floor O(nk)) account for the difference. The
T-010 *kernel* optimisation (cache reorder + C++ candidates) is unchanged; this is
not a regression in the kernel. `dist` input skips the symmetry scan.

## Area 6 — MaxMean (RLTS tabu loop) — AFTER T-012 (monotonicity scan + branchless P-update)

Throughput (tabu iterations/s) is the metric — MaxMean runs for the full
`maxSeconds`, so more iterations = better solution. Signed n=500 Type-I instance,
`useRL = FALSE` (isolates the tabu loop), 3 s budget, `-O2`.

| case | OLD iters/s | NEW iters/s | speedup |
|------|------------:|------------:|--------:|
| MaxMean, n=500, signed, useRL=FALSE | 786k | 1.16M | 1.47× |

Inner-loop self-time split (VTune, n=500): scan ~69%, P-array update ~31%.
`best_flip`/`best_delta` identical (monotonicity); p-arrays bit-identical
(branchless). All tests green; covr 100% on the new code (157/157 C++, 42/42 R).
