# Driver for focus area #6 — MaxMean RLTS tabu inner loop.
# MaxMean is time-budgeted, so the performance metric is THROUGHPUT: iterations
# of the O(n) tabu inner loop achieved per second. A faster inner loop yields
# more iterations in the budget and hence better solutions.
#
# useRL = FALSE isolates the tabu loop (random construction is O(n) and cheap;
# RL construction is O(n^2) per restart but amortized over a long restart).
#
# bare: ~3 s (maxSeconds = 3) on 2026-06-16
#
# Override the library with lib.loc for a -g profiling build, e.g.:
#   Rscript -e '.libPaths("dev/profiling/.vtune-lib-XXatamp"); source(...)'
libloc <- Sys.getenv("MAXMIN_LIB", unset = NA)
if (is.na(libloc)) {
  suppressMessages(library(MaxMin))
} else {
  suppressMessages(library(MaxMin, lib.loc = libloc))
}

set.seed(5813)
n <- 500L
# Signed dissimilarities (the paper's Type-I benchmark range) so the optimum is
# a proper subset and the tabu loop genuinely explores add/remove moves.
m <- matrix(runif(n * n, -10, 10), n)
d <- (m + t(m)) / 2
diag(d) <- 0

budget <- 3
t0 <- proc.time()[["elapsed"]]
res <- MaxMean(d, maxSeconds = budget, maxIter = Inf, useRL = FALSE)
el <- proc.time()[["elapsed"]] - t0

iters <- attr(res, "iters")
cat(sprintf("n=%d  elapsed=%.2fs  iters=%.0f  throughput=%.3g iters/s  |S|=%d  f=%.4f\n",
            n, el, iters, iters / el, length(res), attr(res, "score")))
