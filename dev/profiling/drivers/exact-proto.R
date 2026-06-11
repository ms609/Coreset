# Prototype + verification harness for ExactMaxMin optimisation (area 4).
#
# Three candidate levers, isolated:
#   L1 cheap-edges : precompute upper-tri (row,col,dist) ONCE; each probe is a
#                    vector threshold, not two fresh n x n logical matrices.
#   L2 sparse-A    : build the packing matrix as a sparse dgCMatrix (2 nnz/edge)
#                    instead of a dense nEdge x n matrix (the n>=683 memory wall).
#   L3 feasibility : reformulate each probe "is alpha(G) >= m?" as a constant-
#                    objective feasibility IP with a cardinality cut sum(x)>=m,
#                    so highs stops at the FIRST feasible point on yes-probes
#                    instead of proving the maximum independent set.
#
# Correctness contract: identical PROVEN optimum to the installed ExactMaxMin.
# (The witness subset may differ; the objective lambda* must not.)

set.seed(5813)
suppressMessages({ library(MaxMin); library(Matrix); library(highs) })
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")

# ---- candidate reimplementation -------------------------------------------
Exact_v2 <- function(d, m, maxSeconds = 600, feasibility = TRUE, sparse = TRUE) {
  t0 <- proc.time()[[3L]]
  d <- as.matrix(d); n <- nrow(d); m <- as.integer(m)
  Elapsed <- function() proc.time()[[3L]] - t0

  ut  <- which(upper.tri(d))                 # L1: precompute ONCE
  rc  <- arrayInd(ut, dim(d))
  ui  <- rc[, 1L]; uj <- rc[, 2L]; ud <- d[ut]
  cand <- sort(unique(ud)); nCand <- length(cand)
  card <- sparseMatrix(i = rep(1L, n), j = seq_len(n), x = 1, dims = c(1L, n))

  verdict <- function(lambda, remaining) {
    if (remaining <= 0) return(list(v = "inconclusive", w = integer(0)))
    e <- ud < lambda; nEdge <- sum(e)
    if (nEdge == 0L) return(list(v = "feasible", w = seq_len(n)))
    ei <- ui[e]; ej <- uj[e]
    if (sparse) {
      A <- sparseMatrix(i = rep(seq_len(nEdge), 2L), j = c(ei, ej),
                        x = 1, dims = c(nEdge, n))
    } else {
      A <- matrix(0, nEdge, n)
      A[cbind(seq_len(nEdge), ei)] <- 1; A[cbind(seq_len(nEdge), ej)] <- 1
    }
    if (feasibility) {
      res <- highs_solve(L = rep.int(0, n), lower = rep.int(0, n),
        upper = rep.int(1, n), A = rbind(A, card),
        lhs = c(rep.int(-Inf, nEdge), m), rhs = c(rep.int(1, nEdge), Inf),
        types = rep.int("I", n), maximum = FALSE,
        control = list(threads = 1L, time_limit = remaining))
      sel <- which(res$primal_solution > 0.5)
      if (length(sel) >= m) {
        sub <- d[sel, sel, drop = FALSE]
        if (!any(sub[upper.tri(sub)] < lambda)) return(list(v = "feasible", w = sel))
      }
      st <- res$status_message
      if (identical(st, "Infeasible")) return(list(v = "infeasible", w = integer(0)))
      return(list(v = "inconclusive", w = integer(0)))
    } else {
      res <- highs_solve(L = rep.int(1, n), lower = rep.int(0, n),
        upper = rep.int(1, n), A = A, lhs = rep.int(-Inf, nEdge),
        rhs = rep.int(1, nEdge), types = rep.int("I", n), maximum = TRUE,
        control = list(threads = 1L, time_limit = remaining))
      sel <- which(res$primal_solution > 0.5)
      ok <- if (length(sel) < 2L) TRUE else {
        sub <- d[sel, sel]; !any(sub[upper.tri(sub)] < lambda) }
      opt <- identical(res$status_message, "Optimal")
      if (ok && length(sel) >= m) return(list(v = "feasible", w = sel))
      if (opt && ok && length(sel) < m) return(list(v = "infeasible", w = integer(0)))
      list(v = "inconclusive", w = integer(0))
    }
  }

  lo <- 2L; hi <- nCand; bestIdx <- 1L; bestW <- seq_len(n); inconcl <- FALSE
  while (lo <= hi) {
    rem <- maxSeconds - Elapsed(); if (rem <= 0) { inconcl <- TRUE; break }
    mid <- (lo + hi) %/% 2L
    vv <- verdict(cand[mid], rem)
    if (vv$v == "feasible") { bestIdx <- mid; bestW <- vv$w; lo <- mid + 1L }
    else if (vv$v == "infeasible") { hi <- mid - 1L }
    else { inconcl <- TRUE; break }
  }
  idx <- sort(bestW[seq_len(m)]); sub <- d[idx, idx]; diag(sub) <- Inf
  list(indices = idx, objective = min(sub), proven = !inconcl,
       time_s = Elapsed(), n = n, m = as.integer(m))
}

# ---- compare on the two smallest target cases ------------------------------
bench <- function(case, k = 10L) {
  pts <- as.matrix(cases[[case]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); n <- nrow(d)
  t <- proc.time()[[3L]]; r0 <- MaxMin::ExactMaxMin(d, k, maxSeconds = 600, progress = FALSE); t0 <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]; rA <- Exact_v2(d, k, feasibility = FALSE, sparse = TRUE); tA <- proc.time()[[3L]] - t
  t <- proc.time()[[3L]]; rB <- Exact_v2(d, k, feasibility = TRUE,  sparse = TRUE); tB <- proc.time()[[3L]] - t
  cat(sprintf("%-16s n=%4d k=%d | installed %.2fs (obj %.5f) | v2-maxIS %.2fs (obj %.5f, %.2fx) | v2-feas %.2fs (obj %.5f, %.2fx) | match=%s\n",
    case, n, k, t0, r0$objective, tA, rA$objective, t0/tA, tB, rB$objective, t0/tB,
    isTRUE(all.equal(r0$objective, rA$objective)) && isTRUE(all.equal(r0$objective, rB$objective)) &&
      r0$proven && rA$proven && rB$proven))
}
for (cs in c("tc22_penguins", "tc11_ionosphere")) bench(cs)
