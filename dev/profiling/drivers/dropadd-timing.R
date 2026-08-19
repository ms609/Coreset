# DropAdd kernel timing cells for interleaved A/B minima (this box spikes
# +/-20-35% on sub-second cells; regress against interleaved minima only).
# Most cells call the kernels DIRECTLY, at a fixed maxIter with an unreachable
# plateau, so every build walks the identical deterministic trajectory
# (objective printed as an identity probe). A kernel row is NOT the cost of a
# DropAdd() call: the wrapper's O(n^2) .AsDistMatrix scan and its O(n) seed are
# both outside it, and the `wrapper` and `END-TO-END` cells below price them.
# farfirst-timing.R times FarFirst() through the public API, so only the
# END-TO-END rows here compare with that area's.
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

# Accepts a kernel result (list) or a public-API selection (attributed
# vector); a cell returning neither reports NA and is timed for wall only.
cell <- function(label, th, reps = 3L) {
  ts <- numeric(3L)
  for (r in 1:3) {
    t0 <- proc.time()[[3L]]
    for (q in seq_len(reps)) o <- th()
    ts[r] <- (proc.time()[[3L]] - t0) / reps
  }
  ob <- if (is.list(o)) o$objective else attr(o, "score")
  it <- if (is.list(o)) o$iters     else attr(o, "iters")
  ob <- if (is.null(ob)) NA_real_    else as.numeric(ob)
  it <- if (is.null(it)) NA_integer_ else as.integer(it)
  cat(sprintf("%-44s best=%8.4fs obj=%.10g iters=%d\n", label, min(ts), ob, it))
  data.frame(cell = label, best = min(ts), obj = ob, iters = it)
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
  # What the kernel-direct cells above EXCLUDE. A caller pays these once per
  # DropAdd() call on a fresh matrix; the intake scan is O(n^2), so on the
  # small-m cells it dwarfs the kernel and a kernel row is not a call cost.
  cell("wrapper intake n=4e3 (.AsDistMatrix)", function()
    Coreset:::.AsDistMatrix(d4), 4L),
  cell("wrapper seed n=4e3 (.PickPoint peripheral)", function()
    Coreset:::.PickPoint(d4, "peripheral"), 20L),
  # Whole public-API calls. farfirst-timing.R times FarFirst() end-to-end, so
  # ONLY these rows are comparable with that area's cells.
  cell("END-TO-END DropAdd(20, d4, plateau=200)", function()
    DropAdd(20L, d4, plateau = 200L), 3L),
  cell("END-TO-END DropAdd(600, d4, plateau=200)", function()
    DropAdd(600L, d4, plateau = 200L), 3L),
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
