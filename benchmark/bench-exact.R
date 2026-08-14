source("benchmark/_init.R")

# Exact node-packing optimum on a small instance (it is NP-hard, so the sizes
# stay tiny). Guards the edge enumeration, the clique search that decides each
# probe, and the binary-search orchestration against regression.
d40 <- BenchDist(40L, 5L, seed = 1L)
Benchmark(ExactMaxMin(8L, d40))

# A size at which the whole solve is the single infeasibility proof at the
# certifying threshold, rather than the warm start and the setup scans.
d400 <- BenchDist(400L, 8L, seed = 3L)
Benchmark(ExactMaxMin(10L, d400))
