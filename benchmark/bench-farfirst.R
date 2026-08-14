source("benchmark/_init.R")

# Gonzalez farthest-first across its three input paths, at sizes representative
# of the application benchmark (hundreds-thousands of candidates, tens chosen).

# Dense-matrix path: ensemble default (best-of-three random-furthest starts) and
# a single named seed strategy.
d2000 <- BenchDist(2000L, 10L, seed = 1L)
Benchmark(FarFirst(20L, d2000))                       # default: O(N) ensemble
Benchmark(FarFirst(20L, d2000, strategy = "diameter"))

# Coordinate (matrix-free) path: never materialises the N x N matrix.
pts5000 <- BenchPoints(5000L, 10L, seed = 2L)
Benchmark(FarFirst(k = 30L, points = pts5000))

# Distance-column oracle path
colFn <- function(i) sqrt(rowSums(sweep(pts5000, 2L, pts5000[i, ]) ^ 2))
Benchmark(FarFirst(30L, colFn, N = 5000L))
