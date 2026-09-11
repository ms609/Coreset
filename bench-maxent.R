source("benchmark/_init.R")

d500 <- BenchDist(500L, 8L, seed = 1L)
Benchmark(MaxEntropy(20L, d500))
Benchmark(MaxEntropy(250L, d500))
