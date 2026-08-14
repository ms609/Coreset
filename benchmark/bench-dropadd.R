source("benchmark/_init.R")

# Dense-matrix path.
d500 <- BenchDist(500L, 8L, seed = 1L)
Benchmark(DropAdd(20L, d500, plateau = 2000L))

# Matrix-free coordinate path, at a size past a comfortable dense matrix.
pts4000 <- BenchPoints(4000L, 8L, seed = 2L)
Benchmark(DropAdd(20L, points = pts4000, plateau = 1000L))
