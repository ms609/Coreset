# Verify the NEW (installed) ExactMaxMin: (1) proven optima bit-match the OLD
# baseline across every ground-truth case x k; (2) the speedup; (3) the RNG
# contract -- objective seed-INDEPENDENT, selection reproducible under set.seed.
suppressMessages(library(MaxMin))
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
old <- readRDS("dev/profiling/oldopt.rds")

rows <- list()
for (i in seq_len(nrow(old))) {
  cs <- old$case[i]; k <- old$k[i]
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts))
  set.seed(1)
  t <- proc.time()[[3L]]
  r <- MaxMin::ExactMaxMin(d, k, maxSeconds = 600)
  el <- proc.time()[[3L]] - t
  rows[[i]] <- data.frame(case = cs, n = nrow(d), k = k,
    old = old$objective[i], new = r$objective,
    match = isTRUE(all.equal(old$objective[i], r$objective)) && r$proven,
    old_s = old$time_s[i], new_s = el, speedup = old$time_s[i] / el,
    stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)
cat("================ NEW vs OLD (all ground-truth cases x k) ================\n")
print(R[, c("case","n","k","old","new","match","old_s","new_s","speedup")],
      row.names = FALSE, digits = 6)
cat(sprintf("\nobjective match (proven & equal): %d / %d\n", sum(R$match), nrow(R)))
cat(sprintf("OLD total %.1fs  ->  NEW total %.1fs   (%.1fx overall)\n",
    sum(R$old_s), sum(R$new_s), sum(R$old_s) / sum(R$new_s)))
cat(sprintf("per-case speedup: min %.1fx  median %.1fx  max %.1fx\n",
    min(R$speedup), median(R$speedup), max(R$speedup)))

cat("\n================ RNG contract ================\n")
pts <- as.matrix(cases[["tc9_iris"]][["points"]]); storage.mode(pts) <- "double"
d <- as.matrix(stats::dist(pts))
set.seed(1); a <- MaxMin::ExactMaxMin(d, 6L)
set.seed(2); b <- MaxMin::ExactMaxMin(d, 6L)
set.seed(1); a2 <- MaxMin::ExactMaxMin(d, 6L)
cat(sprintf("objective seed-independent (seed1 == seed2): %s\n",
    isTRUE(all.equal(a$objective, b$objective))))
cat(sprintf("selection reproducible under same seed (seed1 == seed1): %s\n",
    identical(a$indices, a2$indices)))
