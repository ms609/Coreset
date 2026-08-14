source("benchmark/_init.R")

# Dense-matrix path, either side of m = n/8, where the recompute branch
# switches which triangle it reads. Below it the kernel reads a row of each
# selected column; above it, down a single column.
d500 <- BenchDist(500L, 8L, seed = 1L)
Benchmark(DropAdd(20L, d500, plateau = 2000L))
Benchmark(DropAdd(250L, d500, plateau = 2000L))

# Matrix-free coordinate path, at a size past a comfortable dense matrix.
pts4000 <- BenchPoints(4000L, 8L, seed = 2L)
Benchmark(DropAdd(20L, points = pts4000, plateau = 1000L))
