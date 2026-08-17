# Old-vs-new bit-identity battery for FarFirst (grasp-battery.R pattern).
# Captures selections + scores for every seeding strategy on both the matrix
# and coordinate paths, the default random_furthest ensemble, and a
# column-oracle case; also asserts points≡matrix cross-path identity within
# the build under test. Every lever in this area must keep all of it
# bit-identical.
# Usage: Rscript farfirst-battery.R <out.rds>            # capture
#        Rscript farfirst-battery.R <new.rds> <base.rds> # capture + compare
library(Coreset)
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1L]]
# Optional: run the whole battery at a fixed thread count (results must be
# identical to a serial capture — the round-7 kernels are nCores-invariant).
nc <- Sys.getenv("MAXMIN_BATTERY_CORES", "")
if (nzchar(nc)) options(mc.cores = as.integer(nc))

rec <- function(x) list(idx = as.integer(x), score = as.numeric(attr(x, "score")))

res <- list()
crossBad <- character(0)
matStrategies <- c("peripheral", "diameter", "anti_medoid", "medoid",
                   "rowsum", "rownorm")
for (ds in 1:3) {
  for (N in c(120L, 300L, 1500L)) {
    for (dim in c(2L, 10L)) {
      set.seed(ds)
      pts <- matrix(rnorm(N * dim), ncol = dim)
      d <- as.matrix(dist(pts))
      ks <- if (N <= 300L) c(2L, 8L, N %/% 2L, N - 1L) else c(2L, 8L, 750L)
      strategies <- if (N <= 300L) matStrategies else "peripheral"
      for (k in ks) {
        key <- paste("ff", ds, N, dim, k, sep = "_")
        # Explicit first index (bare pass), both paths.
        pm <- rec(FarFirst(k, d, strategy = 5L))
        pp <- rec(FarFirst(k, points = pts, strategy = 5L))
        res[[paste0(key, "_first_mx")]] <- pm
        res[[paste0(key, "_first_pt")]] <- pp
        if (!identical(pm, pp)) crossBad <- c(crossBad, paste0(key, "_first"))
        for (st in strategies) {
          sm <- rec(FarFirst(k, d, strategy = st))
          sp <- rec(FarFirst(k, points = pts, strategy = st))
          res[[paste0(key, "_", st, "_mx")]] <- sm
          res[[paste0(key, "_", st, "_pt")]] <- sp
          if (!identical(sm, sp)) crossBad <- c(crossBad, paste0(key, "_", st))
        }
        # anti_centroid is coordinate-only.
        if (N <= 300L) {
          res[[paste0(key, "_anti_centroid_pt")]] <-
            rec(FarFirst(k, points = pts, strategy = "anti_centroid"))
        }
        # Default ensemble (RNG: seed set immediately before each call).
        set.seed(100L + ds)
        em <- FarFirst(k, d)
        set.seed(100L + ds)
        ep <- FarFirst(k, points = pts)
        res[[paste0(key, "_ens_mx")]] <- c(rec(em),
          list(win = attr(em, "winning_strategy")))
        res[[paste0(key, "_ens_pt")]] <- c(rec(ep),
          list(win = attr(ep, "winning_strategy")))
        if (!identical(res[[paste0(key, "_ens_mx")]],
                       res[[paste0(key, "_ens_pt")]])) {
          crossBad <- c(crossBad, paste0(key, "_ens"))
        }
      }
    }
  }
}
# Column-oracle path (pure R; guarded here against accidental drift).
set.seed(4)
octs <- matrix(rnorm(80L * 3L), ncol = 3L)
oFn <- function(i) sqrt(colSums((t(octs) - octs[i, ])^2))
res[["oracle_80_4"]] <- rec(FarFirst(4L, oFn, N = 80L, strategy = 1L))

saveRDS(res, outFile)
cat("captured", length(res), "cases ->", outFile, "\n")
cat("cross-path identity :", if (length(crossBad)) "BROKEN" else "OK",
    length(crossBad), "\n")
if (length(crossBad)) cat(paste0("  ", head(crossBad, 10), collapse = "\n"), "\n")

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  stopifnot(identical(sort(names(base)), sort(names(res))))
  bad <- names(res)[!vapply(names(res), function(nm)
    identical(base[[nm]], res[[nm]]), logical(1L))]
  cat("\n=== FARFIRST BATTERY ===\n")
  cat("cases compared      :", length(res), "\n")
  cat("mismatches vs base  :", length(bad), "\n")
  if (length(bad)) cat(paste0("  ", head(bad, 20), collapse = "\n"), "\n")
  cat(if (!length(bad) && !length(crossBad))
    "RESULT: BIT-IDENTICAL\n" else "RESULT: MISMATCH\n")
}
