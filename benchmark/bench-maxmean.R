source("benchmark/_init.R")

set.seed(1)
n <- 80
x <- matrix(runif(n * n, -5, 5), n)
d <- (x + t(x)) / 2          # symmetric, signed
Benchmark(MaxMean(d))
