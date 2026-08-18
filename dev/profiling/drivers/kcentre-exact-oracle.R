# Probe-level oracle for the ExactKCentre cover search: every threshold the
# bisection visits is decided BOTH by `CoverDecide_cpp` and by the
# minimum-cardinality set-cover IP it replaces (`highs`), and the two verdicts
# must agree. Because the search is exhaustive, this is exact agreement, not
# merely verdict-preservation: an "infeasible" from the kernel asserts that no
# cover of size <= k exists, which is exactly what the IP proves when its
# certified minimum cover exceeds k.
#
# Feasible witnesses are validated independently of both: the returned centres
# must cover every point within the probed radius.
#
# Counts only -- no wall-clock claim is made from this script.
#
# Usage: Rscript kcentre-exact-oracle.R [maxN]
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
options(Coreset.progress = FALSE)
stopifnot(requireNamespace("highs", quietly = TRUE),
          requireNamespace("Matrix", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
maxN <- if (length(args) >= 1L) as.integer(args[[1L]]) else 400L

load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
CASES <- c("tc20_zoo", "tc7_ring", "tc9_iris", "tc8_highdim_gaussians",
           "tc3_multiscale", "tc1_uniform", "tc2_two_unequal",
           "tc6_density_gradient", "tc16_sonar", "tc5_clusters_outliers",
           "tc10_glass", "tc4_hierarchical", "tc22_penguins",
           "tc11_ionosphere")
KS <- c(2L, 4L, 6L, 10L)

Cands <- getFromNamespace(".KCentreCandidates", "Coreset")
Decide <- getFromNamespace("CoverDecide_cpp", "Coreset")
AsD   <- getFromNamespace(".AsDistMatrix", "Coreset")

# The set-cover IP this round removes, kept here as the oracle.
IpVerdict <- function(d, n, r, k) {
  cover <- which(d <= r, arr.ind = TRUE)
  A <- Matrix::sparseMatrix(i = cover[, 2L], j = cover[, 1L], x = 1,
                            dims = c(n, n))
  res <- highs::highs_solve(
    L = rep.int(1, n), lower = rep.int(0, n), upper = rep.int(1, n), A = A,
    lhs = rep.int(1, n), rhs = rep.int(Inf, n), types = rep.int("I", n),
    maximum = FALSE, control = list(threads = 1L, time_limit = 300))
  sel <- which(res$primal_solution > 0.5)
  ok <- length(sel) >= 1L && all(apply(d[sel, , drop = FALSE] <= r, 2L, any))
  if (ok && length(sel) <= k) return("feasible")
  if (identical(res$status_message, "Optimal") && ok) return("infeasible")
  "inconclusive"
}

nProbe <- 0L; nDisagree <- 0L; nBadWitness <- 0L
totNodes <- 0; maxNodes <- 0
for (cs in CASES) {
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- AsD(stats::dist(pts)); n <- nrow(d)
  if (n > maxN) next
  cand <- Cands(d)
  for (k in KS) {
    set.seed(1000L + k)
    ws <- KCentre(k, d)
    hi <- findInterval(attr(ws, "radius"), cand); lo <- 1L
    cellNodes <- 0
    while (lo < hi) {
      mid <- (lo + hi) %/% 2L; r <- cand[[mid]]
      got <- Decide(d, r, k, 300)
      ip <- IpVerdict(d, n, r, k)
      nProbe <- nProbe + 1L
      cellNodes <- cellNodes + got[["nodes"]]
      totNodes <- totNodes + got[["nodes"]]
      maxNodes <- max(maxNodes, got[["nodes"]])
      if (!identical(got[["status"]], ip)) {
        nDisagree <- nDisagree + 1L
        cat(sprintf("MISMATCH %s k=%d idx=%d: kernel=%s ip=%s\n",
                    cs, k, mid, got[["status"]], ip))
      }
      if (identical(got[["status"]], "feasible")) {
        w <- got[["witness"]]
        if (length(w) < 1L || length(w) > k || is.unsorted(w, strictly = TRUE) ||
            KCentreRadius(d, w) > r) {
          nBadWitness <- nBadWitness + 1L
          cat(sprintf("BAD WITNESS %s k=%d idx=%d\n", cs, k, mid))
        }
      }
      if (identical(ip, "feasible")) hi <- mid else lo <- mid + 1L
    }
    cat(sprintf("%-22s n=%4d k=%2d  nodes=%7.0f\n", cs, n, k, cellNodes))
  }
}
cat(sprintf("\nprobes %d | verdict mismatches %d | bad witnesses %d\n",
            nProbe, nDisagree, nBadWitness))
cat(sprintf("search nodes: total %.0f, worst probe %.0f\n", totNodes, maxNodes))
