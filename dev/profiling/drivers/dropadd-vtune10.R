# VTune driver, round 10: serial hotspot distribution of the DropAdd kernels
# across both m regimes and both paths, kernel-direct (wrapper validation
# excluded — that cost belongs to PR #5's AllFinite lever).
# Usage: Rscript dropadd-vtune10.R <libdir>
args <- commandArgs(trailingOnly = TRUE)
library(Coreset, lib.loc = args[[1L]])
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

set.seed(5813)
n4 <- 4000L
pts4 <- matrix(rnorm(n4 * 10L), ncol = 10L)
d4 <- as.matrix(dist(pts4))
big2  <- matrix(rnorm(20000L * 2L),  ncol = 2L)
big10 <- matrix(rnorm(20000L * 10L), ncol = 10L)
PL <- .Machine$integer.max

t0 <- proc.time()[[3L]]
for (r in 1:20) o1 <- Coreset:::DropAdd_cpp(d4, 10L, Inf, 1500L, PL, FALSE, -1L)
for (r in 1:8)  o2 <- Coreset:::DropAdd_cpp(d4, 2000L, Inf, 1500L, PL, FALSE, -1L)
for (r in 1:4)  o3 <- Coreset:::DropAdd_points_cpp(big2, 10L, Inf, 1000L, PL, FALSE, -1L)
for (r in 1:3)  o4 <- Coreset:::DropAdd_points_cpp(big10, 10L, Inf, 600L, PL, FALSE, -1L)
o5 <- Coreset:::DropAdd_points_cpp(big10, 10000L, Inf, 0L, PL, FALSE, -1L)
cat(sprintf("Elapsed: %.2f s; objs %.6g %.6g %.6g %.6g %.6g\n",
            proc.time()[[3L]] - t0, o1$objective, o2$objective,
            o3$objective, o4$objective, o5$objective))
