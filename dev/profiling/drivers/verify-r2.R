# verify-r2.R <libpath> <tag>
# Correctness + timing battery for round-2 changes (T-004 FarFirst column reorder,
# T-005 DropAdd seed cache-fix, T-006 Grasp incremental swap). All three are
# designed to be BIT-IDENTICAL to the prior build, so correctness == "results
# match the baseline lib exactly". Writes dev/profiling/verify-r2-<tag>.rds.
args <- commandArgs(trailingOnly = TRUE)
libpath <- normalizePath(args[[1]]); tag <- args[[2]]
suppressMessages({ library(Coreset, lib.loc = libpath); library(bench) })
ms <- function(th, it = 11L)
  as.numeric(bench::mark(th(), iterations = it, check = FALSE,
                         filter_gc = FALSE)$median) * 1000

corr <- list(); timing <- list()

## ---- correctness: bit-identity vs baseline ----
# FarFirst (points + matrix), DropAdd (matrix + points), Grasp.
for (ds in 1:4) {
  for (N in c(120L, 300L)) {
    for (dim in c(2L, 5L)) {
      set.seed(ds); pts <- matrix(rnorm(N * dim), ncol = dim)
      d <- as.matrix(dist(pts))
      for (n in c(2L, 8L, as.integer(N / 2), N - 1L)) {
        k <- paste("ff", ds, N, dim, n, sep = "_")
        corr[[paste0(k, "_pt")]] <- as.integer(FarFirst(m = n, points = pts, method = 1L))
        corr[[paste0(k, "_mx")]] <- as.integer(FarFirst(d, n, method = 1L))
      }
      for (m in c(3L, 10L, as.integer(N / 2))) {
        da_m <- DropAdd(d, m, plateau = 200L)
        da_p <- DropAdd(points = pts, m = m, plateau = 200L)
        corr[[paste("da_mx", ds, N, dim, m, sep = "_")]] <-
          list(idx = as.integer(da_m), s = as.numeric(attr(da_m, "score")))
        corr[[paste("da_pt", ds, N, dim, m, sep = "_")]] <-
          list(idx = as.integer(da_p), s = as.numeric(attr(da_p, "score")))
      }
      for (m in c(5L, 20L, as.integer(N / 3))) {
        for (es in c(3L, 5L)) for (al in c(0.5, 0.8)) {
          set.seed(7)
          g <- Grasp(d, m, plateau = 10L, eliteSize = es, alpha = al)
          corr[[paste("gp", ds, N, dim, m, es, al, sep = "_")]] <-
            list(idx = as.integer(g), z = as.numeric(attr(g, "score")))
        }
      }
    }
  }
}

## ---- timing: the three speedups ----
set.seed(11)
# T-005 DropAdd matrix seed (construction-only via maxIter=0), large n small m.
for (N in c(4000L, 6000L)) {
  pts <- matrix(rnorm(N * 2L), ncol = 2L); d <- as.matrix(dist(pts))
  timing[[paste0("da_seed_n", N)]] <- list(
    what = "DropAdd matrix construction (maxIter=0)", N = N,
    ms = ms(function() DropAdd(d, 10L, maxIter = 0L, plateau = 1L)))
}
# T-006 Grasp large m.
set.seed(12)
for (cfg in list(c(200,50), c(200,100))) {
  N <- cfg[1]; m <- cfg[2]
  pts <- matrix(rnorm(N * 2L), ncol = 2L); d <- as.matrix(dist(pts))
  timing[[paste0("grasp_n", N, "_m", m)]] <- list(
    what = "Grasp", N = N, m = m,
    ms = ms(function() { set.seed(1); Grasp(d, m, plateau = 15L, eliteSize = 5L) }, it = 5L))
}
# T-004 FarFirst points single pass (the column reorder), large N, both dims.
set.seed(13)
for (dim in c(2L, 10L)) {
  N <- 6000L; pts <- matrix(rnorm(N * dim), ncol = dim)
  for (n in c(300L, 1500L, 3000L)) {
    timing[[paste0("ff_pts_d", dim, "_n", n)]] <- list(
      what = "FarFirst points single pass", dim = dim, n = n,
      ms = ms(function() FarFirst(m = n, points = pts, method = 1L)))
  }
}

saveRDS(list(tag = tag, corr = corr, timing = timing),
        sprintf("dev/profiling/verify-r2-%s.rds", tag))
cat(sprintf("[%s] %d correctness cases, %d timing cases\n",
            tag, length(corr), length(timing)))
