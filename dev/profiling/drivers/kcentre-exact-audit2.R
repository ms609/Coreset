# ExactKCentre audit 2: where probe cost concentrates, and how far standard
# set-cover reductions shrink the model. Cost figures here are LOCAL SHARES
# for triage only -- no wall-clock claim is made from this box.
#
# Usage: Rscript kcentre-exact-audit2.R out.rds
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(Coreset, lib.loc = lib) else library(Coreset)
options(Coreset.progress = FALSE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "kcentre-audit2.rds"

load(Sys.getenv("FP_CASES",
                "C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda"))
CELLS <- list(c("tc22_penguins", 10), c("tc11_ionosphere", 10),
              c("tc10_glass", 4), c("tc9_iris", 6), c("tc1_uniform", 10))

MinCover <- getFromNamespace(".MinCoverVerdict", "Coreset")
Cands    <- getFromNamespace(".KCentreCandidates", "Coreset")
AsD      <- getFromNamespace(".AsDistMatrix", "Coreset")

# Standard set-cover reductions on the coverage incidence C (C[i, j] = centre i
# covers point j), applied to fixpoint:
#   point j is redundant   if some other point j' has cover(j') subset of cover(j)
#   centre i is redundant  if some other centre i' covers a superset of centre i
# Reports the surviving model size.
Reduce2 <- function(C) {
  pts <- seq_len(ncol(C)); ctr <- seq_len(nrow(C))
  repeat {
    S <- C[ctr, pts, drop = FALSE]
    cs <- colSums(S)
    M <- crossprod(S + 0)                       # M[j, j'] = |cover(j) & cover(j')|
    sub <- M == rep(cs, each = length(cs))      # cover(col) subset of cover(row)
    diag(sub) <- FALSE
    # keep one representative of each equivalence class: drop j if a strictly
    # earlier-indexed j' has a subset cover
    dropPt <- apply(sub & upper.tri(sub, diag = FALSE) | (sub & !t(sub)), 1L, any)
    pts2 <- pts[!dropPt]
    S <- C[ctr, pts2, drop = FALSE]
    rs <- rowSums(S)
    R <- tcrossprod(S + 0)
    sup <- R == rs                              # covers(row) subset of covers(col)
    diag(sup) <- FALSE
    dropCt <- apply(sup & upper.tri(sup, diag = FALSE) | (sup & !t(sup)), 1L, any)
    ctr2 <- ctr[!dropCt]
    if (length(pts2) == length(pts) && length(ctr2) == length(ctr)) break
    pts <- pts2; ctr <- ctr2
    if (!length(pts) || !length(ctr)) break
  }
  c(nPt = length(pts), nCtr = length(ctr))
}

rows <- list()
for (cl in CELLS) {
  cs <- cl[[1]]; k <- as.integer(cl[[2]])
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- AsD(stats::dist(pts)); n <- nrow(d)
  cand <- Cands(d)
  set.seed(1000L + k)
  ws <- KCentre(k, d)
  hi <- findInterval(attr(ws, "radius"), cand); lo <- 1L; p <- 0L
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L; r <- cand[[mid]]
    t0 <- proc.time()[[3L]]
    v <- MinCover(d, n, r, k, 120)
    ip <- proc.time()[[3L]] - t0
    red <- Reduce2(d <= r)
    p <- p + 1L
    rows[[length(rows) + 1L]] <- data.frame(
      case = cs, n = n, k = k, probe = p, idx = mid, verdict = v$verdict,
      ip_s = ip, redPt = red[["nPt"]], redCtr = red[["nCtr"]],
      stringsAsFactors = FALSE)
    cat(sprintf("%-16s k=%2d probe %2d idx=%6d %-11s ip=%7.3f  red %4d x %4d (of %d)\n",
                cs, k, p, mid, v$verdict, ip, red[["nPt"]], red[["nCtr"]], n))
    if (identical(v$verdict, "feasible")) hi <- mid else lo <- mid + 1L
  }
}
res <- do.call(rbind, rows)
saveRDS(res, out)
cat("\n--- probe cost concentration (share of each cell's probe total) ---\n")
for (cs in unique(res$case)) {
  z <- res[res$case == cs, ]
  cat(sprintf("%-16s total %.2fs; top probe %.0f%%; infeasible share %.0f%%\n",
              cs, sum(z$ip_s), 100 * max(z$ip_s) / sum(z$ip_s),
              100 * sum(z$ip_s[z$verdict == "infeasible"]) / sum(z$ip_s)))
}
