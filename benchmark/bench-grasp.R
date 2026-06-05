source("benchmark/_init.R")

# GRASP with path relinking (compiled kernel). Seeded so every repeat does the
# identical stochastic work, and bounded by max_no_improve so the time is the
# compute, not a budget.
d300 <- BenchDist(300L, 6L, seed = 1L)
Benchmark(GraspPR(d300, 15L, max_no_improve = 100L, elite_size = 8L, seed = 1L))
