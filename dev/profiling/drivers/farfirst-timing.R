# FarFirst kernel timing cells (single pass, strategy = explicit first), for
# interleaved A/B minima (this box spikes +/-20-35% on sub-second cells).
# Usage: Rscript farfirst-timing.R <out.rds> [nThreads]
library(Coreset)
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1L]]
if (length(args) >= 2L) options(mc.cores = as.integer(args[[2L]]))
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

set.seed(5813)
N <- 6000L
pts2  <- matrix(rnorm(N * 2L),  ncol = 2L)
pts10 <- matrix(rnorm(N * 10L), ncol = 10L)
d10   <- as.matrix(dist(pts10))
big2  <- matrix(rnorm(100000L * 2L),  ncol = 2L)
big10 <- matrix(rnorm(100000L * 10L), ncol = 10L)

cell <- function(label, th, reps) {
  ts <- numeric(3L)
  for (r in 1:3) {
    t0 <- proc.time()[[3L]]
    for (q in seq_len(reps)) o <- th()
    ts[r] <- (proc.time()[[3L]] - t0) / reps
  }
  cat(sprintf("%-34s best=%8.4fs score=%.6g\n", label, min(ts),
              attr(o, "score")))
  data.frame(cell = label, best = min(ts), score = as.numeric(attr(o, "score")))
}

tab <- rbind(
  cell("matrix N=6e3 dim=10 k=3000", function()
    FarFirst(3000L, d10, strategy = 5L), 4L),
  cell("points N=6e3 dim=2 k=3000", function()
    FarFirst(3000L, points = pts2, strategy = 5L), 6L),
  cell("points N=6e3 dim=10 k=3000", function()
    FarFirst(3000L, points = pts10, strategy = 5L), 3L),
  cell("points N=1e5 dim=2 k=1000", function()
    FarFirst(1000L, points = big2, strategy = 5L), 3L),
  cell("points N=1e5 dim=10 k=1000", function()
    FarFirst(1000L, points = big10, strategy = 5L), 2L),
  cell("ens-anchors N=6e3 dim=10 k=300", function()
    FarFirst(300L, points = pts10,
             strategy = c("diameter", "anti_medoid", "rownorm")), 1L),
  cell("ens-anchors-matrix N=6e3 k=300", function()
    FarFirst(300L, d10,
             strategy = c("diameter", "anti_medoid", "rownorm")), 1L),
  # Default ensemble: nSeeds restarts from random-furthest seeds. The seeds are
  # O(N) each, so the restarts are essentially the whole call — the shape the
  # restart-parallel kernels address. Seeded per call: the draw must not vary
  # between the arms of an A/B.
  cell("ens-default N=6e3 dim=10 k=3000", function() {
    set.seed(99); FarFirst(3000L, d10)
  }, 1L),
  cell("ens-default-pts N=6e3 d10 k=3000", function() {
    set.seed(99); FarFirst(3000L, points = pts10)
  }, 1L),
  cell("ens-default nSeeds=8 N=6e3 k=3000", function() {
    set.seed(99); FarFirst(3000L, d10, nSeeds = 8L)
  }, 1L)
)
saveRDS(tab, out)
