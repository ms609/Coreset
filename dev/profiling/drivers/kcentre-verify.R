# OLD-vs-NEW verification for the CDSh optimisation (T-010).
# Usage: Rscript kcentre-verify.R <lib.loc> <tag>
# Saves dev/profiling/.kc-<tag>.rds = list(res by k, ms at k=20).
args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1]]; tag <- args[[2]]
suppressMessages(library(MaxMin, lib.loc = lib))
set.seed(5813)
g <- 12L; per <- 167L
ctr <- matrix(rnorm(g * 10L, sd = 6), ncol = 10L)
pts <- ctr[rep(seq_len(g), each = per), ] + matrix(rnorm(g * per * 10L), ncol = 10L)
d <- as.matrix(dist(pts)); n <- nrow(d)

res <- list()
for (k in c(5L, 20L, 50L)) {
  r <- KCentre(d, k)
  res[[as.character(k)]] <- list(idx = as.integer(r), radius = attr(r, "radius"))
}
t0 <- proc.time()
for (i in 1:6) rr <- KCentre(d, 20L)
ms <- 1000 * (proc.time() - t0)[["elapsed"]] / 6

saveRDS(list(res = res, ms = ms), sprintf("dev/profiling/.kc-%s.rds", tag))
cat(sprintf("%s: per-call=%.1fms  k20 radius=%.4f nidx=%d\n",
            tag, ms, res[["20"]]$radius, length(res[["20"]]$idx)))
