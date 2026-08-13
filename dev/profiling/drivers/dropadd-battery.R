# Old-vs-new trajectory-identity battery for DropAdd (grasp-battery.R
# pattern). Calls the kernels DIRECTLY with want_trace = TRUE so every case
# records the full drop/add sequence alongside indices, objective, secondary
# and iters — any tie-break slip in a fused argmax shows up as a trajectory
# divergence even when the final selection happens to coincide. Axes:
# matrix path (Euclidean, TIE-DENSE integer-valued, and asymmetric inputs;
# default and explicit seeds; small and large m), points path (dim 1/2/7/10),
# plus two within-build invariants that need no baseline:
#   * cross-path: matrix vs points kernels with the same explicit seed on
#     d = dist(pts) must produce the identical trajectory bit-for-bit;
#   * oracle twin: the pure-R .DropAddFromColumn must reproduce the matrix
#     kernel's result exactly (the R twin is this round's frozen semantic
#     reference).
# Usage: Rscript dropadd-battery.R <out.rds>            # capture
#        Rscript dropadd-battery.R <new.rds> <base.rds> # capture + compare
library(MaxMin)
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1L]]
cat("lib:", dirname(system.file(package = "MaxMin")), "\n")

recM <- function(dmat, m, plateau, seed0, maxIter = 100000L) {
  o <- MaxMin:::DropAdd_cpp(dmat, m, Inf, maxIter, plateau, TRUE, seed0)
  list(idx = as.integer(o$indices), obj = as.numeric(o$objective),
       sec = as.numeric(o$secondary), iters = as.integer(o$iters),
       drops = as.integer(o$drops), adds = as.integer(o$adds))
}
recP <- function(pts, m, plateau, seed0, maxIter = 100000L) {
  o <- MaxMin:::DropAdd_points_cpp(pts, m, Inf, maxIter, plateau, TRUE, seed0)
  list(idx = as.integer(o$indices), obj = as.numeric(o$objective),
       sec = as.numeric(o$secondary), iters = as.integer(o$iters),
       drops = as.integer(o$drops), adds = as.integer(o$adds))
}

res <- list()
invBad <- character(0)

for (ds in 1:3) {
  for (n in c(60L, 300L, 800L)) {
    set.seed(ds * 1000L + n)
    pts <- matrix(rnorm(n * 3L), ncol = 3L)
    dEu  <- as.matrix(dist(pts))
    dTie <- matrix(as.double(sample.int(5L, n * n, TRUE)), n, n)
    dTie <- pmax(dTie, t(dTie))          # symmetric, heavily tied
    diag(dTie) <- 0
    dAsy <- matrix(runif(n * n), n, n)   # asymmetric (tolerated input)
    diag(dAsy) <- 0
    ms <- unique(c(2L, 5L, n %/% 4L, n %/% 2L, n - 1L))
    for (m in ms) {
      key <- paste("da", ds, n, m, sep = "_")
      res[[paste0(key, "_eu")]]  <- recM(dEu,  m, 30L, -1L)
      res[[paste0(key, "_tie")]] <- recM(dTie, m, 30L, -1L)
      res[[paste0(key, "_asy")]] <- recM(dAsy, m, 30L, -1L)
      res[[paste0(key, "_eu_s3")]] <- recM(dEu, m, 30L, 2L)  # explicit seed
      # Points kernel (its own anti-centroid default seed) + cross-path
      # trajectory identity at an explicit shared seed.
      res[[paste0(key, "_pt")]]    <- recP(pts, m, 30L, -1L)
      pm <- recM(dEu, m, 30L, 4L)
      pp <- recP(pts, m, 30L, 4L)
      res[[paste0(key, "_pt_s5")]] <- pp
      if (!identical(pm, pp)) invBad <- c(invBad, paste0(key, "_cross"))
    }
    # One deep-plateau cell per shape (long trajectory).
    res[[paste0("da_", ds, "_", n, "_deep")]] <-
      recM(dEu, max(2L, n %/% 5L), 2000L, -1L)
  }
}

# Higher-dim points cells (dim block shapes if fills get blocked).
for (dim in c(1L, 2L, 7L, 10L)) {
  set.seed(77L + dim)
  pts <- matrix(rnorm(400L * dim), ncol = dim)
  for (m in c(2L, 40L, 200L)) {
    res[[paste0("da_dim", dim, "_", m)]] <- recP(pts, m, 30L, -1L)
    dm <- as.matrix(dist(pts))
    pm <- recM(dm, m, 30L, 7L)
    pp <- recP(pts, m, 30L, 7L)
    if (!identical(pm, pp)) {
      invBad <- c(invBad, paste0("da_dim", dim, "_", m, "_cross"))
    }
  }
}

# Oracle twin: the pure-R path must reproduce the matrix kernel exactly.
for (ds in 1:2) {
  set.seed(300L + ds)
  n <- 120L
  pts <- matrix(rnorm(n * 3L), ncol = 3L)
  dmat <- as.matrix(dist(pts))
  colFn <- function(i) dmat[, i]
  for (m in c(3L, 25L)) {
    ro <- MaxMin:::.DropAddFromColumn(colFn, n, m, first = 6L,
                                      plateau = 30L, trace = TRUE)
    co <- MaxMin:::DropAdd_cpp(dmat, m, Inf, 100000L, 30L, TRUE, 5L)
    twin <- identical(sort(as.integer(ro$indices)),
                      sort(as.integer(co$indices))) &&
      identical(as.numeric(ro$objective), as.numeric(co$objective)) &&
      identical(as.numeric(ro$secondary), as.numeric(co$secondary)) &&
      identical(as.integer(ro$iters), as.integer(co$iters)) &&
      identical(as.integer(ro$drops), as.integer(co$drops)) &&
      identical(as.integer(ro$adds), as.integer(co$adds))
    res[[paste0("da_oracle_", ds, "_", m)]] <-
      list(idx = as.integer(ro$indices), obj = as.numeric(ro$objective),
           sec = as.numeric(ro$secondary), iters = as.integer(ro$iters),
           drops = as.integer(ro$drops), adds = as.integer(ro$adds))
    if (!twin) invBad <- c(invBad, paste0("da_oracle_", ds, "_", m))
  }
}

saveRDS(res, outFile)
cat("captured", length(res), "cases ->", outFile, "\n")
cat("within-build invariants:",
    if (length(invBad)) "BROKEN" else "OK", length(invBad), "\n")
if (length(invBad)) cat(paste0("  ", head(invBad, 10), collapse = "\n"), "\n")

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  stopifnot(identical(sort(names(base)), sort(names(res))))
  bad <- names(res)[!vapply(names(res), function(nm)
    identical(base[[nm]], res[[nm]]), logical(1L))]
  cat("\n=== DROPADD BATTERY ===\n")
  cat("cases compared      :", length(res), "\n")
  cat("mismatches vs base  :", length(bad), "\n")
  if (length(bad)) cat(paste0("  ", head(bad, 20), collapse = "\n"), "\n")
  cat(if (!length(bad) && !length(invBad))
    "RESULT: BIT-IDENTICAL\n" else "RESULT: MISMATCH\n")
}
