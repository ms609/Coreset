source("benchmark/_init.R")

# GRASP with path relinking (compiled kernel). set.seed() before the run for a
# reproducible trajectory, and bounded by `plateau` so the time is the compute,
# not a budget.
d300 <- BenchDist(300L, 6L, seed = 1L)
set.seed(1)
Benchmark(Grasp(d300, 15L, plateau = 100L, eliteSize = 8L))
