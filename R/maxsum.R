# maxsum.R
#
# Solver for the Max-Sum Diversity Problem (MDP, "maximum diversity"): given an
# n x n symmetric distance matrix and a target size k, find the k-subset S
# maximising the TOTAL pairwise distance T = sum_{i < j in S} d(i, j). This is
# the max-SUM counterpart of the max-MIN p-dispersion solved by ExactMaxMin().
#
# Method: exact MILP via the Kuo, Glover & Dhir (1993) per-node linearisation,
# the most efficient integer program for the MDP \insertCite{Kuo1993}{Coreset}.
# The quadratic objective sum_{i<j} d_ij x_i x_j = (1/2) sum_i x_i a_i, with
# a_i = sum_j d_ij x_j the contribution of i, is linearised by one continuous
# variable w_i = x_i a_i per point:
#   maximise (1/2) sum_i w_i
#   s.t.  sum_i x_i = k
#         w_i <= U_i x_i            (w_i = 0 unless i is selected)
#         w_i <= sum_j d_ij x_j     (w_i <= a_i)
#         x in {0,1}, w >= 0
# where U_i = sum of the k largest d_ij is an upper bound on a_i. Under
# maximisation the two upper bounds force w_i = a_i exactly when i is selected,
# so no lower-bound constraints are needed. The program has only 2n variables
# (n binary + n continuous) and 1 + 2n constraints -- far smaller and tighter
# than the O(n^2) edge linearisation, which makes it practical to n ~ 100.
#
# The search is floored by a fast multi-start 1-swap local search (the standard
# MDP neighbourhood with an incremental contribution array), which both warms
# the lower bound and is returned when the MILP cannot prove optimality within
# the budget. As the MDP is NP-hard this is, like ExactMaxMin(), a small-
# instance tool.
#
# Solver: the `highs` package (HiGHS MILP backend; CRAN).

# ----- heuristic core -------------------------------------------------------

# Total pairwise distance of a selection -- the max-sum objective.
.MaxSumScore <- function(d, idx) {
  s <- d[idx, idx, drop = FALSE]
  sum(s[upper.tri(s)])
}

# Greedy C2/GMA construction: seed at the most peripheral point (largest row
# sum), then repeatedly add the point maximising the summed distance to the
# selection. Returns a sorted k-subset.
.MaxSumGreedy <- function(d, k, seed) {
  n <- nrow(d)
  chosen <- integer(k)
  chosen[1L] <- seed
  contrib <- d[, seed]
  for (t in seq_len(k - 1L) + 1L) {
    contrib[chosen[seq_len(t - 1L)]] <- -Inf
    chosen[t] <- which.max(contrib)
    contrib <- contrib + d[, chosen[t]]
  }
  sort(chosen)
}

# 1-swap (interchange) local search to a local optimum, using an incremental
# contribution array gain[i] = sum_{j in S} d[i, j]. A swap (u in S, v not in S)
# changes the objective by gain[v] - gain[u] - d[u, v]; the best improving swap
# is applied and gain updated in O(n). Returns the improved sorted subset.
.MaxSumLocalSearch <- function(d, k, idx) {
  n <- nrow(d)
  inS <- logical(n)
  inS[idx] <- TRUE
  gain <- colSums(d[idx, , drop = FALSE])
  repeat {
    sel <- which(inS)
    uns <- which(!inS)
    # Best improving interchange over all (u in S) x (v not in S).
    delta <- outer(gain[uns], gain[sel], "-") - d[uns, sel, drop = FALSE]
    best <- which.max(delta)
    if (delta[best] <= 1e-12) {
      break
    }
    vi <- uns[(best - 1L) %% length(uns) + 1L]
    ui <- sel[(best - 1L) %/% length(uns) + 1L]
    inS[ui] <- FALSE
    inS[vi] <- TRUE
    gain <- gain + d[, vi] - d[, ui]
  }
  sort(which(inS))
}

# Multi-start greedy + local search: best subset over `nStart` peripheral and
# random seeds. Deterministic given the session RNG. Returns a sorted k-subset.
.MaxSumHeuristic <- function(d, k, nStart = 8L) {
  n <- nrow(d)
  rowSum <- rowSums(d)
  seeds <- unique(c(which.max(rowSum),
                    order(rowSum, decreasing = TRUE)[seq_len(min(nStart, n))],
                    sample.int(n, min(nStart, n))))
  best <- NULL
  bestVal <- -Inf
  for (s in seeds) {
    idx <- .MaxSumLocalSearch(d, k, .MaxSumGreedy(d, k, s))
    v <- .MaxSumScore(d, idx)
    if (v > bestVal) {
      bestVal <- v
      best <- idx
    }
  }
  best
}

# ----- exact MILP -----------------------------------------------------------

#' Exact Maximum Diversity Problem (max-sum) solution
#'
#' `ExactMaxSum()` finds the optimal solution to the Max-Sum Diversity Problem
#' (the "maximum diversity problem"): select the `k`-subset of points
#' maximising the **total** pairwise distance it contains. It is the max-sum
#' counterpart of [ExactMaxMin()] (which maximises the *minimum* pairwise
#' distance), solved by per-node integer-program linearisation
#' \insertCite{Kuo1993}{Coreset}. As the problem is NP-hard it is feasible only
#' for small sets.
#'
#' The optimum is floored by a multi-start 1-swap local search, which warms the
#' lower bound and is returned when the MILP cannot prove optimality within
#' `maxSeconds` -- so the result is always at least a strong heuristic incumbent.
#'
#' @inheritParams ExactMaxMin
#' @return `ExactMaxSum()` returns an integer vector of length `k` (sorted
#'   ascending) with class `"MaxSumSelection"`, carrying attributes:
#'   \describe{
#'     \item{score}{Achieved total pairwise distance within the selection. When
#'       `proven` is `TRUE` this is the optimum; otherwise a lower bound.}
#'     \item{proven}{Logical: `TRUE` if optimality was certified within the
#'       budget.}
#'     \item{time_s, N, k}{Wall-clock seconds, instance size, target size.}
#'   }
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(20), ncol = 2)
#' ExactMaxSum(3L, dist(pts))
#' @export
ExactMaxSum <- function(k, d, maxSeconds = 60, warmStart = NULL) {
  t0 <- proc.time()[[3L]]
  if (!requireNamespace("highs", quietly = TRUE)) { # nocov start
    stop("The `highs` package is required for ExactMaxSum(). ",
         "Install it with install.packages(\"highs\").")
  } # nocov end

  d <- .ExactAsMatrix(d)
  n <- nrow(d)
  k <- as.integer(k)
  if (is.na(k) || k < 2L || k > n) {
    stop("`k` must satisfy 2 <= k <= nrow(d)")
  }
  Elapsed <- function() proc.time()[[3L]] - t0

  # Heuristic floor (also a valid incumbent if the MILP times out). An optional
  # caller `warmStart` joins the comparison.
  lsIdx <- .MaxSumHeuristic(d, k)
  lsVal <- .MaxSumScore(d, lsIdx)
  if (!is.null(warmStart)) {
    ws <- tryCatch(sort(unique(as.integer(warmStart))), error = function(e) integer(0))
    if (length(ws) == k && ws[1L] >= 1L && ws[k] <= n) {
      wv <- .MaxSumScore(d, ws)
      if (wv > lsVal) { lsVal <- wv; lsIdx <- ws }
    }
  }

  # Kuo--Glover--Dhir MILP. Variables x_1..x_n (binary), w_1..w_n (continuous);
  # U_i bounds a_i above by the sum of the k largest distances from i.
  U <- vapply(seq_len(n), function(i) {
    sum(sort(d[i, -i], decreasing = TRUE)[seq_len(min(k, n - 1L))])
  }, double(1))
  ri <- rep(1L, n); ci <- seq_len(n); vi <- rep(1, n)            # sum x = k
  rA <- 1L + seq_len(n)                                          # w_i <= U_i x_i
  ri <- c(ri, rA, rA); ci <- c(ci, n + seq_len(n), seq_len(n)); vi <- c(vi, rep(1, n), -U)
  rB <- 1L + n + seq_len(n)                                      # w_i <= sum_j d_ij x_j
  ri <- c(ri, rB, rep(rB, each = n))
  ci <- c(ci, n + seq_len(n), rep(seq_len(n), times = n))
  vi <- c(vi, rep(1, n), -as.vector(d))
  A <- Matrix::sparseMatrix(i = ri, j = ci, x = vi, dims = c(1L + 2L * n, 2L * n))

  res <- highs::highs_solve(
    L       = c(rep(0, n), rep(0.5, n)),
    lower   = rep(0, 2L * n),
    upper   = c(rep(1, n), U),
    A       = A,
    lhs     = c(k, rep(-Inf, 2L * n)),
    rhs     = c(k, rep(0, 2L * n)),
    types   = c(rep("I", n), rep("C", n)),
    maximum = TRUE,
    control = list(threads = 1L, time_limit = max(0, maxSeconds - Elapsed()))
  )

  milpIdx <- sort(which(res$primal_solution[seq_len(n)] > 0.5))
  proven <- identical(res$status_message, "Optimal") && length(milpIdx) == k
  milpVal <- if (length(milpIdx) == k) .MaxSumScore(d, milpIdx) else -Inf

  # Return the better of the certified/incumbent MILP solution and the local
  # search. Optimality is claimed only when the MILP proved it.
  if (milpVal >= lsVal) {
    idx <- milpIdx; score <- milpVal
  } else {
    idx <- lsIdx; score <- lsVal; proven <- FALSE
  }

  # Return:
  structure(idx, score = score, proven = proven, time_s = Elapsed(),
            N = n, k = as.integer(k), producer = "ExactMaxSum",
            class = "MaxSumSelection")
}
