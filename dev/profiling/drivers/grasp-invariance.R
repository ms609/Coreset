# T-021 determinism gate: a given seed must return a bit-identical result at
# EVERY thread count (nCores trades wall-clock only). Random draws happen
# solely on the main thread in fixed-size batches and batches merge in
# iteration order, so any divergence here is a bug, not noise.
library(MaxMin)
cat("lib:", dirname(system.file(package = "MaxMin")), "\n")

shapes <- list(
  list(n = 2000L, dim = 10L, k = 100L, pl = 64L),
  list(n = 500L,  dim = 2L,  k = 50L,  pl = 64L),
  list(n = 200L,  dim = 10L, k = 10L,  pl = 128L)
)
bad <- 0L
for (sh in shapes) {
  set.seed(11)
  d <- as.matrix(dist(matrix(rnorm(sh$n * sh$dim), sh$n)))
  res <- lapply(c(1L, 2L, 8L), function(nc) {
    set.seed(1)
    o <- MaxMin:::Grasp_cpp(d, sh$k, sh$pl, .Machine$integer.max, 10L, 0.8,
                            Inf, nc)
    list(idx = as.integer(o$indices), obj = as.numeric(o$objective),
         iters = as.numeric(o$iters), pr = as.numeric(o$pr_calls))
  })
  ok <- identical(res[[1L]], res[[2L]]) && identical(res[[1L]], res[[3L]])
  cat(sprintf("n=%d k=%d plateau=%d nCores {1,2,8}: %s\n", sh$n, sh$k, sh$pl,
              if (ok) "INVARIANT" else "*** DIFFERS ***"))
  if (!ok) bad <- bad + 1L
}
cat(if (bad == 0L) "RESULT: NCORES-INVARIANT\n" else "RESULT: MISMATCH\n")
