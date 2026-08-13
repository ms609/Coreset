# VTune driver, round 9: serial hotspot distribution of the FarFirst kernels
# AND the O(N^2) anchor primitives (round 7's driver covered the passes only).
# Single Gonzalez passes (strategy = explicit first) at the design-regime
# sizes, plus the serial anchors ensemble whose RowSums/RowSqSums/Diameter
# primitives round 7 parallelised but never line-profiled.
# Usage: Rscript farfirst-vtune9.R <libdir>
args <- commandArgs(trailingOnly = TRUE)
library(MaxMin, lib.loc = args[[1L]])
cat("lib:", dirname(system.file(package = "MaxMin")), "\n")

set.seed(5813)
N <- 6000L
pts10 <- matrix(rnorm(N * 10L), ncol = 10L)
d10   <- as.matrix(dist(pts10))
big2  <- matrix(rnorm(100000L * 2L),  ncol = 2L)
big10 <- matrix(rnorm(100000L * 10L), ncol = 10L)

t0 <- proc.time()[[3L]]
for (r in 1:8) {
  o1 <- FarFirst(3000L, d10, strategy = 5L)
}
for (r in 1:6) {
  o2 <- FarFirst(1000L, points = big2,  strategy = 5L)
}
for (r in 1:3) {
  o3 <- FarFirst(1000L, points = big10, strategy = 5L)
}
for (r in 1:3) {
  o4 <- FarFirst(300L, points = pts10,
                 strategy = c("diameter", "anti_medoid", "rownorm"))
}
for (r in 1:3) {
  o5 <- FarFirst(300L, d10,
                 strategy = c("diameter", "anti_medoid", "rownorm"))
}
cat(sprintf("Elapsed: %.2f s; scores %.6g %.6g %.6g %.6g %.6g\n",
            proc.time()[[3L]] - t0, attr(o1, "score"), attr(o2, "score"),
            attr(o3, "score"), attr(o4, "score"), attr(o5, "score")))
