# Driver: MaxEntropy() maxdet selector internal-step breakdown.
# The one O(n^3) cost is the eigendecomposition(s) in the R wrapper; everything
# else is O(n^2) / O(n^2 k). The wrapper does TWO eigens (negMass values-only +
# repair full). We measure the cost and the gain from merging them into one.
# Distance matrix is a synthetic non-Euclidean one: eigen FLOPs depend on SIZE
# not values, and a random symmetric "distance" makes the RBF kernel non-PSD so
# the clip-repair path is genuinely exercised (negMass > 0).
.libPaths(c("C:/Users/pjjg18/GitHub/MaxMin/.agent-lib", .libPaths()))
suppressMessages(library(MaxMin))
set.seed(5813)
tmN <- function(expr, R) { e <- new.env(); median(replicate(3, system.time(for (i in seq_len(R)) eval(expr, e))[["elapsed"]] / R)) }

mergedPrep <- function(K, tol = 1e-9) {
  ks <- (K + t(K)) / 2
  e <- eigen(ks, symmetric = TRUE)
  lam <- e$values; neg <- lam[lam < -tol]
  negMass <- if (length(lam)) sum(abs(neg)) / sum(abs(lam)) else 0
  lam[lam < 0] <- 0
  list(kp = e$vectors %*% (lam * t(e$vectors)), negMass = negMass)
}

for (n in c(1000L, 1500L, 2000L)) {
  M <- matrix(runif(n * n), n); D <- (M + t(M)); diag(D) <- 0     # non-Euclidean
  K <- MaxMin:::.MaxEntropyKernel(D)
  R <- if (n <= 1000) 4L else 2L
  t_neg    <- tmN(quote(MaxMin:::.MaxEntropyNegMass(K)), R)
  t_repair <- tmN(quote(MaxMin:::.MaxEntropyRepair(K, "clip")), R)
  t_merged <- tmN(quote(mergedPrep(K)), R)
  t_full   <- tmN(quote(MaxEntropy(10L, D)), R)
  cur <- t_neg + t_repair
  cat(sprintf("n=%4d | negMass %.3f + repair %.3f = %.3f s  ->  merged %.3f s  | prep gain %4.1f%% | full=%.3f s, est full gain %4.1f%%\n",
              n, t_neg, t_repair, cur, t_merged, 100*(cur - t_merged)/cur, t_full, 100*(cur - t_merged)/t_full))
}
