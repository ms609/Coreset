# Correctness battery for ExactKCentre over the ground-truth grid:
# 14 cases x k in {2, 4, 6, 10} (56 cells). Captures per cell: proven flag,
# radius, witness, wall time; validates the witness independently of whatever
# produced it (length <= k, sorted unique indices in range, achieved covering
# radius bit-equal to the reported `radius`).
#
# Between builds the witness is FREE to differ (ties admit several optima; the
# contract, inherited from ExactMaxMin's round-11 steer, is "any provably
# optimal solution"). The gates are `radius` bit-identity, `proven` identity,
# and independent witness validity. Radius identity is witness-independent: a
# proven optimum's witness achieves exactly the smallest feasible candidate,
# since its own achieved radius is itself a feasible candidate.
#
# Usage:
#   Rscript kcentre-exact-grid.R out.rds            # capture
#   Rscript kcentre-exact-grid.R new.rds base.rds   # capture + compare
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
options(Coreset.progress = FALSE)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)

casePath <- Sys.getenv("FP_CASES",
                       "C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
load(casePath)
grid <- expand.grid(
  case = c("tc20_zoo", "tc7_ring", "tc9_iris", "tc8_highdim_gaussians",
           "tc3_multiscale", "tc1_uniform", "tc2_two_unequal",
           "tc6_density_gradient", "tc16_sonar", "tc5_clusters_outliers",
           "tc10_glass", "tc4_hierarchical", "tc22_penguins",
           "tc11_ionosphere"),
  k = c(2L, 4L, 6L, 10L),
  stringsAsFactors = FALSE
)

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
  r <- ExactKCentre(k, d, maxSeconds = 600)
  wall <- proc.time()[[3L]] - t0
  idx <- as.integer(r)
  radius <- attr(r, "radius")
  valid <- length(idx) >= 1L && length(idx) <= k &&
    !is.unsorted(idx, strictly = TRUE) &&
    idx[[1L]] >= 1L && idx[[length(idx)]] <= n &&
    identical(KCentreRadius(d, idx), radius)
  cells[[i]] <- list(case = cs, n = n, k = k, radius = radius,
                     proven = isTRUE(attr(r, "proven")), valid = valid,
                     idx = idx, wall = wall)
  cat(sprintf("%-22s n=%4d k=%2d  radius=%.6f proven=%s valid=%s  %.2fs\n",
              cs, n, k, radius, cells[[i]]$proven, valid, wall))
}
saveRDS(cells, args[[1L]])
cat(sprintf("\ntotal %.2f s; all valid: %s\n",
            sum(vapply(cells, function(z) z$wall, numeric(1))),
            all(vapply(cells, function(z) z$valid, logical(1)))))

if (length(args) >= 2L) {
  base <- readRDS(args[[2L]])
  bad <- 0L
  for (i in seq_along(cells)) {
    a <- cells[[i]]; b <- base[[i]]
    if (!identical(a$radius, b$radius) || !identical(a$proven, b$proven) ||
        !a$valid) {
      bad <- bad + 1L
      cat(sprintf("MISMATCH %-22s k=%2d  radius %.17g vs %.17g  proven %s/%s\n",
                  a$case, a$k, a$radius, b$radius, a$proven, b$proven))
    }
  }
  cat(sprintf("\ncompared %d cells; %d mismatches\n", length(cells), bad))
}
