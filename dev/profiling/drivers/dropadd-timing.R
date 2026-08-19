# DropAdd kernel timing cells for interleaved A/B minima (this box spikes
# +/-20-35% on sub-second cells; regress against interleaved minima only).
# Kernels are called DIRECTLY: the R wrapper's .AsDistMatrix validation cost
# belongs to PR #5's AllFinite lever, not to this area. All cells run a
# fixed maxIter with an unreachable plateau, so every build walks the
# identical deterministic trajectory (objective printed as an identity
# probe).
# The matrix cells pass the SAME peripheral seed DropAdd() computes, so they
# time the trajectory a caller actually walks; it is computed once, outside
# the timed thunks. Matrix figures recorded before round 14 timed the kernel's
# O(n^2) max-row-sum seed instead and are NOT comparable with these. The
# points cells pass -1L: the anti-centroid fallback IS that path's production
# default.
# Usage: Rscript dropadd-timing.R <out.rds>
library(Coreset)
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1L]]
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

set.seed(5813)
n4 <- 4000L
pts4 <- matrix(rnorm(n4 * 10L), ncol = 10L)
d4 <- as.matrix(dist(pts4))
big2  <- matrix(rnorm(20000L * 2L),  ncol = 2L)
big10 <- matrix(rnorm(20000L * 10L), ncol = 10L)
PL <- .Machine$integer.max
# Production protocol: DropAdd()'s matrix default seed, 0-based.
s4 <- as.integer(Coreset:::.PickPoint(d4, "peripheral")) - 1L

cell <- function(label, th, reps = 3L) {
  ts <- numeric(3L)
  for (r in 1:3) {
    t0 <- proc.time()[[3L]]
    for (q in seq_len(reps)) o <- th()
    ts[r] <- (proc.time()[[3L]] - t0) / reps
  }
  cat(sprintf("%-36s best=%8.4fs obj=%.10g iters=%d\n", label, min(ts),
              o$objective, o$iters))
  data.frame(cell = label, best = min(ts), obj = o$objective,
             iters = o$iters)
}

tab <- rbind(
  cell("matrix n=4e3 m=10 construct", function()
    Coreset:::DropAdd_cpp(d4, 10L, Inf, 0L, PL, FALSE, s4), 6L),
  cell("matrix n=4e3 m=2000 construct", function()
    Coreset:::DropAdd_cpp(d4, 2000L, Inf, 0L, PL, FALSE, s4), 4L),
  cell("matrix n=4e3 m=10 search1500", function()
    Coreset:::DropAdd_cpp(d4, 10L, Inf, 1500L, PL, FALSE, s4), 4L),
  cell("matrix n=4e3 m=2000 search1500", function()
    Coreset:::DropAdd_cpp(d4, 2000L, Inf, 1500L, PL, FALSE, s4), 3L),
  # Either side of m = n/8, where the recompute branch switches read order.
  cell("matrix n=4e3 m=400 search1500", function()
    Coreset:::DropAdd_cpp(d4, 400L, Inf, 1500L, PL, FALSE, s4), 3L),
  cell("matrix n=4e3 m=600 search1500", function()
    Coreset:::DropAdd_cpp(d4, 600L, Inf, 1500L, PL, FALSE, s4), 3L),
  cell("points n=2e4 d2 m=10 search1000", function()
    Coreset:::DropAdd_points_cpp(big2, 10L, Inf, 1000L, PL, FALSE, -1L), 3L),
  cell("points n=2e4 d10 m=10 search600", function()
    Coreset:::DropAdd_points_cpp(big10, 10L, Inf, 600L, PL, FALSE, -1L), 3L),
  cell("points n=2e4 d10 m=1e4 construct", function()
    Coreset:::DropAdd_points_cpp(big10, 10000L, Inf, 0L, PL, FALSE, -1L), 1L),
  cell("points n=2e4 d10 m=1e4 search600", function()
    Coreset:::DropAdd_points_cpp(big10, 10000L, Inf, 600L, PL, FALSE, -1L), 1L)
)
saveRDS(tab, out)
