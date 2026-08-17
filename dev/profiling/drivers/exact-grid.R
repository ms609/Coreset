# Correctness battery for ExactMaxMin over the ground-truth grid:
# 14 cases x k in {2, 4, 6, 10} (56 cells). Captures per cell: proven flag,
# score, witness, wall time; validates the witness independently of the solver
# (length k, sorted unique indices in range, achieved min pairwise distance
# bit-equal to `score`).
#
# Between builds the witness is FREE to differ (ties admit several optima; the
# user's contract is "any provably optimal solution"): the gates are `score`
# bit-identity, `proven` identity, and independent witness validity.
#
# Usage:
#   Rscript exact-grid.R out.rds            # capture
#   Rscript exact-grid.R new.rds base.rds   # capture + compare
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
options(Coreset.progress = FALSE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)

load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
grid <- expand.grid(
  case = c("tc20_zoo", "tc7_ring", "tc9_iris", "tc8_highdim_gaussians",
           "tc3_multiscale", "tc1_uniform", "tc2_two_unequal",
           "tc6_density_gradient", "tc16_sonar", "tc5_clusters_outliers",
           "tc10_glass", "tc4_hierarchical", "tc22_penguins",
           "tc11_ionosphere"),
  k = c(2L, 4L, 6L, 10L),
  stringsAsFactors = FALSE
)

.Scv <- function(d, idx) {
  s <- d[idx, idx]
  diag(s) <- Inf
  min(s)
}

cells <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  cs <- grid$case[[i]]
  k <- grid$k[[i]]
  pts <- as.matrix(cases[[cs]][["points"]])
  storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts))
  n <- nrow(d)
  set.seed(1000L + i)
  t0 <- proc.time()[[3L]]
  r <- ExactMaxMin(k, d, maxSeconds = 600)
  wall <- proc.time()[[3L]] - t0
  idx <- as.integer(r)
  score <- attr(r, "score")
  valid <- length(idx) == k && !is.unsorted(idx, strictly = TRUE) &&
    idx[[1L]] >= 1L && idx[[k]] <= n && identical(.Scv(d, idx), score)
  cells[[i]] <- list(case = cs, n = n, k = k, score = score,
                     proven = isTRUE(attr(r, "proven")), valid = valid,
                     idx = idx, wall = wall)
  cat(sprintf("%-22s n=%4d k=%2d  score=%.6f proven=%s valid=%s  %.2fs\n",
              cs, n, k, score, cells[[i]]$proven, valid, wall))
}
saveRDS(cells, args[[1L]])

allValid <- all(vapply(cells, `[[`, logical(1), "valid"))
allProven <- all(vapply(cells, `[[`, logical(1), "proven"))
cat(sprintf("\nAll witnesses valid: %s   All proven: %s   Total wall: %.1fs\n",
            allValid, allProven, sum(vapply(cells, `[[`, numeric(1), "wall"))))

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  stopifnot(length(base) == length(cells))
  bad <- 0L
  for (i in seq_along(cells)) {
    a <- cells[[i]]
    b <- base[[i]]
    if (!identical(a$score, b$score) || !identical(a$proven, b$proven)) {
      bad <- bad + 1L
      cat(sprintf("MISMATCH %s k=%d: score %.10g/%.10g proven %s/%s\n",
                  a$case, a$k, a$score, b$score, a$proven, b$proven))
    }
  }
  cat(sprintf("\nCompared %d cells vs %s: %d mismatches\n",
              length(cells), args[[2L]], bad))
  cat(sprintf("Wall: base %.1fs -> new %.1fs\n",
              sum(vapply(base, `[[`, numeric(1), "wall")),
              sum(vapply(cells, `[[`, numeric(1), "wall"))))
}
