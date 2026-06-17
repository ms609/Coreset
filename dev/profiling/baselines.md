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

| case | metric | ms |
|------|--------|---:|
| matrix construction, n=4000, m=10 (maxIter=0) | median | 10.3 |
| matrix construction, n=6000, m=10 (maxIter=0) | median | 22.4 |
| matrix search, n=4000, per-iter | — | ~0.06 |
| points search, n=20000 dim=2, per-iter | — | ~0.23 |

(Pre-fix construction: 116 ms / 293 ms → 11.3× / 13.1×.)

## Area 3 — Grasp — AFTER T-007 (base_z min-edge witness hoist)

| case | metric | ms |
|------|--------|---:|
| n=200, m=50, eliteSize=5, plateau=15 | median | 16.4 |
| n=200, m=100, eliteSize=5, plateau=15 | median | 39.5 |

(T-006→T-007: 22.2→16.4 ms / 85.0→39.5 ms → 1.34× / 2.07×, 416/416 bit-identical.
 Cumulative vs pre-T-006: 315→16.4 ms / 1329→39.5 ms → 19.2× / 33.6×.)

## Area 1 — FarFirst single pass — AFTER T-004 (column reorder, AT-LIMIT)

| case | metric | ms |
|------|--------|---:|
| points, dim=2, N=6000, n=3000 | median | 40.4 |
| points, dim=10, N=6000, n=3000 | median | 108.8 |

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
