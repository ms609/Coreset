source("benchmark/_init.R")

# Euclidean distances in eight dimensions give a numerically positive-definite
# RBF kernel (the Cholesky fast path); d^1.2 breaks the negative-type property
# so the kernel is indefinite and the partial eigendecomposition runs.
d500 <- BenchDist(500L, 8L, seed = 1L)
d500Indefinite <- d500 ^ 1.2

Benchmark(MaxEntropy(20L, d500))
Benchmark(MaxEntropy(250L, d500))
Benchmark(MaxEntropy(20L, d500Indefinite))
Benchmark(MaxEntropy(20L, d500Indefinite, repair = "shift"))
