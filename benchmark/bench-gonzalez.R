source("benchmark/_init.R")

# Gonzalez farthest-first across its three input paths, at sizes representative
# of the application benchmark (hundreds-thousands of candidates, tens chosen).

# Dense-matrix path: ensemble default (best-of-three random-furthest starts) and
# a single named seed strategy.
d2000 <- BenchDist(2000L, 10L, seed = 1L)
Benchmark(Gonzalez(d2000, 20L))                       # default: O(N) ensemble
Benchmark(Gonzalez(d2000, 20L, seed = "diameter"))

# Coordinate (matrix-free) path: never materialises the N x N matrix.
pts5000 <- BenchPoints(5000L, 10L, seed = 2L)
Benchmark(Gonzalez(points = pts5000, n = 30L))

# Distance-column oracle path (as used by TreeSearch on tree sets): pass a
# closure as `d` and Gonzalez() reads one column at a time, never a matrix.
# Here the column is recomputed from coordinates on the fly.
colFn <- function(i) sqrt(rowSums(sweep(pts5000, 2L, pts5000[i, ]) ^ 2))
Benchmark(Gonzalez(colFn, 30L, N = 5000L))
