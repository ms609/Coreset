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

## Area 3 — GraspPR — AFTER T-006 (incremental swap)

| case | metric | ms |
|------|--------|---:|
| n=200, m=50, eliteSize=5, plateau=15 | median | 22.2 |
| n=200, m=100, eliteSize=5, plateau=15 | median | 85.0 |

(Pre-fix: 315 ms / 1329 ms → 14.2× / 15.6×.)

## Area 1 — FarFirst single pass — AFTER T-004 (column reorder, AT-LIMIT)

| case | metric | ms |
|------|--------|---:|
| points, dim=2, N=6000, n=3000 | median | 40.4 |
| points, dim=10, N=6000, n=3000 | median | 108.8 |
