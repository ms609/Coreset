source("benchmark/_init.R")

# GRASP with path relinking (compiled kernel).
d300 <- BenchDist(300L, 6L, seed = 1L)
set.seed(1)
Benchmark(Grasp(15L, d300, plateau = 100L, eliteSize = 8L))
