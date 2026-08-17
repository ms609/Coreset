# Cores-vs-speedup curve for Grasp (T-021). Results are nCores-invariant
# (asserted by grasp-invariance.R), so the objective column must not move.
# Local relative figures for choosing between implementations; real scaling
# curves belong on Hamilton.
library(Coreset)
cat("lib:", dirname(system.file(package = "Coreset")), "\n")
set.seed(11)
d <- as.matrix(dist(matrix(rnorm(2000 * 10), 2000)))
REPS <- 3L
for (pl in c(8L, 256L)) {
  for (nc in c(1L, 2L, 4L, 8L)) {
    ts <- numeric(REPS)
    for (r in seq_len(REPS)) {
      set.seed(1)
      t0 <- proc.time()[[3L]]
      o <- Coreset:::Grasp_cpp(d, 100L, pl, .Machine$integer.max, 10L, 0.8,
                              Inf, nc)
      ts[r] <- proc.time()[[3L]] - t0
    }
    cat(sprintf("plateau=%-3d nCores=%d best=%6.3fs med=%6.3fs obj=%.6f iters=%d\n",
                pl, nc, min(ts), median(ts), o$objective,
                as.integer(o$iters)))
  }
}
