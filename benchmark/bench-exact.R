source("benchmark/_init.R")

# Exact node-packing optimum on a small instance (it is NP-hard, so the sizes
# stay tiny). The MILP solve itself is in `highs`; this guards our edge
# enumeration, IP assembly and binary-search orchestration against regression.
if (requireNamespace("highs", quietly = TRUE)) {
  d40 <- BenchDist(40L, 5L, seed = 1L)
  Benchmark(ExactMaxMin(8L, d40))
}
