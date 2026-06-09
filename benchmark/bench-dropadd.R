source("benchmark/_init.R")

# DropAdd tabu search, both paths, under a fixed deterministic iteration budget
# (max_no_improve) so the measured time reflects the work, not a wall clock.

# Dense-matrix path.
d500 <- BenchDist(500L, 8L, seed = 1L)
Benchmark(DropAdd(d500, 20L, max_no_improve = 2000L))

# Matrix-free coordinate path, at a size past a comfortable dense matrix.
pts4000 <- BenchPoints(4000L, 8L, seed = 2L)
Benchmark(DropAdd(points = pts4000, m = 20L, plateau = 1000L))
