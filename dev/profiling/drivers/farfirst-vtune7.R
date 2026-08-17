# VTune driver, round 7: serial hotspot distribution of the FarFirst kernels.
# Single Gonzalez passes (strategy = explicit first) so the kernels dominate:
# matrix path at N=6000, points path at dim {2, 10}, plus one large-N points
# case (the path's design regime). Usage: Rscript farfirst-vtune7.R <libdir>
# bare: ~4 s of kernel work after the one-off dist() setup, 2026-08-13
args <- commandArgs(trailingOnly = TRUE)
library(Coreset, lib.loc = args[[1L]])
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

set.seed(5813)
N <- 6000L
pts2  <- matrix(rnorm(N * 2L),  ncol = 2L)
pts10 <- matrix(rnorm(N * 10L), ncol = 10L)
d10   <- as.matrix(dist(pts10))
big   <- matrix(rnorm(60000L * 2L), ncol = 2L)

t0 <- proc.time()[[3L]]
for (r in 1:12) {
  o1 <- FarFirst(3000L, d10, strategy = 5L)
  o2 <- FarFirst(3000L, points = pts2,  strategy = 5L)
  o3 <- FarFirst(3000L, points = pts10, strategy = 5L)
}
for (r in 1:6) {
  o4 <- FarFirst(1000L, points = big, strategy = 5L)
}
cat(sprintf("Elapsed: %.2f s; scores %.6g %.6g %.6g %.6g\n",
            proc.time()[[3L]] - t0, attr(o1, "score"), attr(o2, "score"),
            attr(o3, "score"), attr(o4, "score")))
