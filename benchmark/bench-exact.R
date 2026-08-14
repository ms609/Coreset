source("benchmark/_init.R")

# Exact node-packing optimum on a small instance (it is NP-hard, so the sizes
# stay tiny). Guards the edge enumeration, the clique search that decides each
# probe, and the binary-search orchestration against regression.
d40 <- BenchDist(80L, 5L, seed = 1L)
Benchmark(ExactMaxMin(8L, d40))

# A size at which the whole solve is the single infeasibility proof at the
# certifying threshold, rather than the warm start and the setup scans.
d1200 <- BenchDist(1200L, 12L, seed = 3L)
Benchmark(ExactMaxMin(12L, d1200))
