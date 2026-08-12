# Old-vs-new bit-identity battery for Grasp (T-007 `oldopt.rds` pattern).
# Usage: Rscript battery.R <out.rds>            # capture
#        Rscript battery.R <new.rds> <base.rds> # capture + compare
library(MaxMin)
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[[1L]]

grid <- expand.grid(
  n         = c(60L, 200L, 500L),
  dim       = c(2L, 10L),
  k         = c(3L, 10L, 50L),
  eliteSize = c(1L, 4L, 10L),
  alpha     = c(0, 0.8, 1),
  plateau   = c(1L, 8L, 40L),
  seed      = 1:3,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
grid <- grid[grid$k < grid$n, ]

res <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  set.seed(1000L + g$seed)
  d <- as.matrix(dist(matrix(rnorm(g$n * g$dim), g$n)))
  set.seed(g$seed)
  o <- MaxMin:::Grasp_cpp(d, g$k, g$plateau, .Machine$integer.max,
                          g$eliteSize, as.double(g$alpha), Inf)
  res[[i]] <- list(
    idx = as.integer(o$indices), obj = as.numeric(o$objective),
    iters = as.numeric(o$iters), pr = as.numeric(o$pr_calls)
  )
  # Also exercise the R reference on the small cases: re-proves parity directly.
  if (g$n <= 60L && g$k <= 10L) {
    set.seed(g$seed)
    r <- MaxMin:::.Grasp_R(g$k, d, plateau = g$plateau,
                           eliteSize = g$eliteSize, alpha = as.double(g$alpha))
    res[[i]]$rIdx <- as.integer(r)
    res[[i]]$rObj <- as.numeric(attr(r, "score"))
  }
}
saveRDS(list(grid = grid, res = res), outFile)
cat("captured", nrow(grid), "cells ->", outFile, "\n")

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  stopifnot(identical(dim(base$grid), dim(grid)))
  bad <- character(0)
  for (i in seq_along(res)) {
    a <- base$res[[i]]; b <- res[[i]]
    lab <- paste0("cell", i, " ", paste(names(grid), unlist(grid[i, ]),
                                       sep = "=", collapse = " "))
    if (!identical(a$idx, b$idx))                bad <- c(bad, paste(lab, "INDICES"))
    if (!identical(a$obj, b$obj))                bad <- c(bad, paste(lab, "OBJECTIVE"))
    if (!identical(a$iters, b$iters))            bad <- c(bad, paste(lab, "ITERS"))
    if (!identical(a$pr, b$pr))                  bad <- c(bad, paste(lab, "PR_CALLS"))
  }
  # R-reference parity, checked within the new build itself.
  rp <- 0L; rbad <- character(0)
  for (i in seq_along(res)) {
    b <- res[[i]]
    if (!is.null(b$rIdx)) {
      rp <- rp + 1L
      if (!identical(b$rIdx, b$idx) || !identical(b$rObj, b$obj)) {
        rbad <- c(rbad, paste0("cell", i))
      }
    }
  }
  cat("\n=== BATTERY ===\n")
  cat("cells compared      :", length(res), "\n")
  cat("mismatches vs base  :", length(bad), "\n")
  if (length(bad)) cat(paste0("  ", head(bad, 20), collapse = "\n"), "\n")
  cat("R-vs-C++ parity     :", rp - length(rbad), "/", rp, "\n")
  if (length(rbad)) cat("  R parity broke:", paste(head(rbad, 20), collapse = ", "), "\n")
  cat(if (!length(bad) && !length(rbad)) "RESULT: BIT-IDENTICAL\n" else "RESULT: MISMATCH\n")
}
