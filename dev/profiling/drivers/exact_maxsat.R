# Clique kernel with MaxSAT branching-set reduction (ThresholdDecide_cpp).
# Real certifying / last-feasible probes of ExactMaxMin cells (vowel, breast
# cancer, vehicle, penguins, pima), built by FurthestPoint data; the fixture
# path comes from FP_PROBES. Pure C++ from the Rcpp boundary inward, so no
# profvis step. The two pima k=100 probes run to a fixed cap: they carry the
# long-propagation-chain regime (up to 99 classes) the others lack.
lib <- Sys.getenv("CORESET_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
Decide <- Coreset:::ThresholdDecide_cpp
probes <- readRDS(Sys.getenv("FP_PROBES", "probes.rds"))
maxsat <- as.logical(Sys.getenv("FP_MAXSAT", "TRUE"))
reps <- as.integer(Sys.getenv("FP_REPS", "8"))
capBig <- as.numeric(Sys.getenv("FP_CAP", "1.5"))
set.seed(5813)
t0 <- proc.time()
tot <- 0
for (p in probes) {
  big <- p$k >= 100
  r <- NULL
  t1 <- proc.time()[["elapsed"]]
  for (i in seq_len(if (big) 1L else reps)) {
    r <- Decide(p$hi, p$hj, p$n, p$k, if (big) capBig else 600, 1L, maxsat)
  }
  el <- proc.time()[["elapsed"]] - t1
  tot <- tot + r$nodes * (if (big) 1 else reps)
  cat(sprintf("%-34s %-12s %8.3f s %12.0f nodes\n", p$label, r$status, el, r$nodes))
}
cat("Elapsed:", round((proc.time() - t0)[["elapsed"]], 2), "s; nodes", tot, "\n")
