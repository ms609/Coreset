library(Coreset)
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1L]]
cat("lib:", dirname(system.file(package = "Coreset")), "\n")

shapes <- list(
  list(n = 2000L, dim = 10L, k = 100L, pl = 8L),
  list(n = 2000L, dim = 10L, k = 100L, pl = 64L),
  list(n = 2000L, dim = 10L, k = 100L, pl = 256L),
  list(n = 2000L, dim = 10L, k = 10L,  pl = 64L),
  list(n = 500L,  dim = 10L, k = 50L,  pl = 64L)
)
REPS <- 3L
tab <- data.frame()
for (sh in shapes) {
  set.seed(11)
  d <- as.matrix(dist(matrix(rnorm(sh$n * sh$dim), sh$n)))
  ts <- numeric(REPS)
  for (r in seq_len(REPS)) {
    set.seed(1)
    t0 <- proc.time()[[3]]
    o <- Coreset:::Grasp_cpp(d, sh$k, sh$pl, .Machine$integer.max, 10L, 0.8, Inf)
    ts[r] <- proc.time()[[3]] - t0
  }
  lab <- sprintf("n=%d k=%d plateau=%d", sh$n, sh$k, sh$pl)
  tab <- rbind(tab, data.frame(shape = lab, best = min(ts), med = median(ts),
                               iters = as.numeric(o$iters),
                               obj = as.numeric(o$objective)))
  cat(sprintf("%-26s best=%7.3fs med=%7.3fs iters=%d\n", lab, min(ts),
              median(ts), as.integer(o$iters)))
}
saveRDS(tab, out)
