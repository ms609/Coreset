# VTune driver, round 6: serial hotspot distribution of the 0.0.0.9008 kernel.
# Canonical coreset shape only (n=2000, dim=10, k=100, plateau=256) so the
# collection is dominated by the path the baselines quote (~0.8 s/call serial).
# Usage: Rscript grasp-vtune6.R <libdir> [reps]
# bare: ~5 s (6 reps) on 2026-08-13
args <- commandArgs(trailingOnly = TRUE)
libdir <- args[[1L]]
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 6L
library(Coreset, lib.loc = libdir)
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

set.seed(11)
d <- as.matrix(dist(matrix(rnorm(2000L * 10L), 2000L)))
t0 <- proc.time()[[3L]]
for (r in seq_len(reps)) {
  set.seed(1)
  o <- Coreset:::Grasp_cpp(d, 100L, 256L, .Machine$integer.max, 10L, 0.8, Inf,
                          1L)
}
cat(sprintf("Elapsed: %.2f s over %d reps; iters=%d obj=%.17g\n",
            proc.time()[[3L]] - t0, reps, as.integer(o$iters),
            as.numeric(o$objective)))
