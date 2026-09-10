# Reduction-potential audit for ExactKCentre (T-011), the covering dual of
# ExactMaxMin. COUNTS ONLY -- no wall-clock claim is made from this script; all
# timing for this area belongs on Hamilton.
#
# For every threshold the shipped bisection actually probes, this records:
#   verdict        -- what the highs covering IP said (feasible / infeasible)
#   greedy         -- size of a greedy set cover (<= k would settle "feasible"
#                     with no IP; the covering analogue of round 11's greedy
#                     k-clique shortcut on the max-min side)
#   confClique     -- size of a clique found in the CONFLICT graph on points
#                     (j ~ j' iff no centre covers both). Points in a conflict
#                     clique need pairwise distinct centres, so a clique of
#                     size >= k+1 proves the threshold infeasible with no IP.
#                     Sound for arbitrary symmetric d -- it does NOT assume the
#                     triangle inequality, unlike the "> 2r apart" shortcut.
#   nnz            -- non-zeros in the covering incidence (IP model size)
#
# Usage: Rscript kcentre-exact-audit.R out.rds [maxN]
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
options(Coreset.progress = FALSE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "kcentre-audit.rds"
maxN <- if (length(args) >= 2L) as.integer(args[[2L]]) else 400L

load(Sys.getenv("FP_CASES",
                "C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda"))
CASES <- c("tc20_zoo", "tc7_ring", "tc9_iris", "tc8_highdim_gaussians",
           "tc3_multiscale", "tc1_uniform", "tc2_two_unequal",
           "tc6_density_gradient", "tc16_sonar", "tc5_clusters_outliers",
           "tc10_glass", "tc4_hierarchical", "tc22_penguins",
           "tc11_ionosphere")
KS <- c(2L, 4L, 6L, 10L)

MinCover <- getFromNamespace(".MinCoverVerdict", "Coreset")
Cands    <- getFromNamespace(".KCentreCandidates", "Coreset")
AsD      <- getFromNamespace(".AsDistMatrix", "Coreset")
Decide   <- getFromNamespace("ThresholdDecide_cpp", "Coreset")

# Greedy minimum set cover over the coverage incidence `C` (C[i, j] = centre i
# covers point j). Deterministic: ties break to the lowest centre index.
GreedyCover <- function(C, k) {
  n <- ncol(C)
  uncov <- rep(TRUE, n)
  picked <- integer(0)
  # Cap the walk: we only care whether greedy lands at or below k.
  while (any(uncov) && length(picked) <= k + 2L) {
    gain <- as.integer(C[, uncov, drop = FALSE] %*% rep(1L, sum(uncov)))
    b <- which.max(gain)
    if (gain[[b]] == 0L) break
    picked <- c(picked, b)
    uncov[uncov] <- !C[b, uncov]
  }
  if (any(uncov)) NA_integer_ else length(picked)
}

# A clique in the conflict graph, via the ExactMaxMin kernel: ask for size
# `k + 1`. Returns the target size when one exists, 0 when the exhaustive
# search proves none does, NA when it timed out.
ConflictClique <- function(C, k, budget = 5) {
  # co[j, j'] = number of centres covering both points; 0 => conflict edge.
  co <- crossprod(C + 0)
  conf <- co == 0
  conf[!upper.tri(conf)] <- FALSE
  e <- which(conf, arr.ind = TRUE)
  if (!nrow(e)) return(0L)
  r <- Decide(as.integer(e[, 1L]), as.integer(e[, 2L]), ncol(C),
              as.integer(k + 1L), budget)
  switch(r$status, feasible = k + 1L, infeasible = 0L, NA_integer_)
}

rows <- list()
for (cs in CASES) {
  pts <- as.matrix(cases[[cs]][["points"]])
  storage.mode(pts) <- "double"
  d <- AsD(stats::dist(pts))
  n <- nrow(d)
  if (n > maxN) { cat(sprintf("%-22s n=%4d  SKIP (> maxN)\n", cs, n)); next }
  cand <- Cands(d)
  for (k in KS) {
    set.seed(1000L + k)
    ws <- KCentre(k, d)
    hi <- findInterval(attr(ws, "radius"), cand)
    wsIdx <- hi
    lo <- 1L
    optIdx <- hi
    nProbe <- 0L
    while (lo < hi) {
      mid <- (lo + hi) %/% 2L
      r <- cand[[mid]]
      C <- d <= r
      v <- MinCover(d, n, r, k, 60)
      g <- GreedyCover(C, k)
      cq <- ConflictClique(C, k)
      nProbe <- nProbe + 1L
      rows[[length(rows) + 1L]] <- data.frame(
        case = cs, n = n, k = k, probe = nProbe, idx = mid,
        nCand = length(cand), wsIdx = wsIdx,
        verdict = v$verdict, greedy = g, confClique = cq,
        nnz = sum(C), stringsAsFactors = FALSE)
      if (identical(v$verdict, "feasible")) { hi <- mid; optIdx <- mid
      } else lo <- mid + 1L
    }
    cat(sprintf("%-22s n=%4d k=%2d  ws_idx=%5d opt_idx=%5d probes=%2d\n",
                cs, n, k, length(cand), wsIdx, optIdx, nProbe))
  }
}
res <- do.call(rbind, rows)
saveRDS(res, out)
cat("\n--- verdict mix ---\n"); print(table(res$verdict))
cat("\n--- would greedy settle the feasible probes? (greedy <= k) ---\n")
f <- res[res$verdict == "feasible", ]
print(table(settled = !is.na(f$greedy) & f$greedy <= f$k))
cat("\n--- would a conflict clique settle the infeasible probes? (>= k+1) ---\n")
i <- res[res$verdict == "infeasible", ]
print(table(settled = !is.na(i$confClique) & i$confClique >= i$k + 1L))
