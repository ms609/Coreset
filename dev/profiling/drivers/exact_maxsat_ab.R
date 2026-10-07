# Interleaved A/B of the clique kernel with and without MaxSAT absorption on the real probe fixture
# (FP_PROBES, as exact_maxsat.R): per probe, each arm's median of FP_ROUNDS
# timed repeats, arms alternating within every round. Probes that do not
# complete within FP_CAP report nodes per second instead.
lib <- Sys.getenv("CORESET_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
Decide <- Coreset:::ThresholdDecide_cpp
probes <- readRDS(Sys.getenv("FP_PROBES", "probes.rds"))
arms <- as.logical(strsplit(Sys.getenv("FP_ARMS", "FALSE,TRUE"), ",")[[1]])  # maxsat
rounds <- as.integer(Sys.getenv("FP_ROUNDS", "5"))
cap <- as.numeric(Sys.getenv("FP_CAP", "2"))
inner <- as.integer(Sys.getenv("FP_INNER", "20"))   # repeats per timing (10 ms clock)
out <- list()
for (p in probes) {
  tm <- matrix(NA_real_, rounds, length(arms)); nd <- tm; st <- character(length(arms))
  for (r in seq_len(rounds)) for (a in seq_along(arms)) {
    t0 <- proc.time()[["elapsed"]]
    big <- p$k >= 100
    for (i in seq_len(if (big) 1L else inner)) {
      res <- Decide(p$hi, p$hj, p$n, p$k, cap, 1L, arms[a])
    }
    tm[r, a] <- (proc.time()[["elapsed"]] - t0) / (if (big) 1L else inner)
    nd[r, a] <- res$nodes; st[a] <- res$status
  }
  done <- all(st != "inconclusive")
  row <- data.frame(probe = p$label, status = paste(unique(st), collapse = "/"))
  for (a in seq_along(arms)) {
    row[[paste0("ms", arms[a])]] <- if (done) median(tm[, a]) else median(nd[, a] / tm[, a])
  }
  row$metric <- if (done) "secs" else "nodes/s"
  for (a in seq_along(arms)) row[[paste0("n", arms[a])]] <- nd[1, a]
  if (length(unique(st)) > 1) stop("verdicts differ on ", p$label)
  out[[length(out) + 1]] <- row
  print(row, row.names = FALSE); flush.console()
}
