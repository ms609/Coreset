# exact.R
#
# Exact solver for the Max-Min Diversity Problem (MMDP = discrete
# p-dispersion): given an n x n symmetric distance matrix and a target size
# m, find the m-subset S maximising T_k = min_{i != j in S} d(i, j).
#
# Method: iterated node-packing (Sayyady & Fathi 2016, EJOR 253(1):216-225).
# The MaxMin optimum is the largest threshold lambda over the achieved
# distinct pairwise distances such that the threshold graph
#   G(lambda) = (V, { (i, j) : d(i, j) < lambda })
# admits an independent set of size >= m. An independent set in G(lambda) is
# a set of points all pairwise >= lambda apart -- i.e. a feasible m-subset
# with min-distance >= lambda. We binary-search the sorted distinct
# distances (descending feasibility) and, at each tested lambda, solve a
# maximum-independent-set feasibility IP
#   maximize sum_i x_i  s.t.  x_i + x_j <= 1 for every edge (i, j);  x in {0,1}
# declaring lambda feasible iff the optimum (the independence number) >= m.
#
# Used in the manuscript ONLY as an external ground-truth reference on small
# instances -- it is NP-hard and not a scalable method.
#
# Solver: the `highs` package (Florian Schwendinger; HiGHS MILP backend;
# CRAN). Chosen over Rglpk/Rsymphony for speed and a self-contained binary.

# ----- helpers --------------------------------------------------------------

# Local coercion to a square numeric matrix. Deliberately does not call the
# package internal `.AsDistMatrix` (which may move); a dist object or square
# numeric matrix is all this solver needs.
.ExactAsMatrix <- function(d) {
  if (inherits(d, "dist")) {
    return(as.matrix(d))
  }
  if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a `dist` object or a square numeric matrix")
  }
  d
}

# Solve one maximum-independent-set IP on the threshold graph G(lambda) and
# classify the result against the target size m.
#
# Builds one packing constraint x_i + x_j <= 1 per edge (i, j) with
# d(i, j) < lambda, then maximises sum x with binary variables. The returned
# witness is validated independently of the solver status: x is rounded, the
# selected set is checked to contain no G(lambda) edge, and its size is
# compared to m. This makes the classification robust to solver
# status-code or time-limit quirks.
#
# Returns list(verdict, witness) where verdict is one of:
#   "feasible"     -- a validated independent set of size >= m was found
#                     (witness = its vertex indices; min-distance >= lambda),
#   "infeasible"   -- the IP was solved to proven optimality and the maximum
#                     independent set has size < m (witness = integer(0)),
#   "inconclusive" -- the budget expired before either could be established
#                     (witness = integer(0)).
.MaxISVerdict <- function(d, n, lambda, m, timeLimit) {
  if (!is.finite(timeLimit) || timeLimit <= 0) { # nocov start
    return(list(verdict = "inconclusive", witness = integer(0)))
  } # nocov end
  # Edges of G(lambda): unordered pairs with d(i, j) < lambda. `<` (strict)
  # is exact here because lambda is itself an achieved pairwise distance and
  # both d and lambda come from the same matrix -- no rounding intervenes.
  edge <- which(d < lambda & upper.tri(d), arr.ind = TRUE)
  nEdge <- nrow(edge)

  if (nEdge == 0L) { # nocov start
    # Empty graph: every vertex is independent, so alpha = n. Feasible
    # whenever m <= n (guaranteed by the caller's guard). No solve needed.
    return(list(verdict = "feasible", witness = seq_len(n)))
  } # nocov end

  # One row per edge: x_i + x_j <= 1.
  A <- matrix(0, nrow = nEdge, ncol = n)
  A[cbind(seq_len(nEdge), edge[, 1L])] <- 1
  A[cbind(seq_len(nEdge), edge[, 2L])] <- 1

  res <- highs::highs_solve(
    L       = rep.int(1, n),
    lower   = rep.int(0, n),
    upper   = rep.int(1, n),
    A       = A,
    lhs     = rep.int(-Inf, nEdge),
    rhs     = rep.int(1, nEdge),
    types   = rep.int("I", n),                # integer var on [0,1] = binary
    maximum = TRUE,
    control = list(
      threads    = 1L,                        # determinism
      time_limit = timeLimit
    )
  )

  # Independent validation of the witness -- never trust the status alone.
  sel <- which(res$primal_solution > 0.5)
  isValidIndependent <- if (length(sel) < 2L) {
    TRUE  # nocov                         # 0 or 1 vertex: trivially independent
  } else {
    sub <- d[sel, sel, drop = FALSE]
    !any(sub[upper.tri(sub)] < lambda)
  }

  if (isValidIndependent && length(sel) >= m) {
    return(list(verdict = "feasible", witness = sel))
  }

  # Infeasibility is provable ONLY when the IP reached optimality and the
  # certified maximum independent set is still smaller than m. A time-limit
  # hit with a too-small (or empty) incumbent proves nothing.
  optimal <- identical(res$status_message, "Optimal")
  if (optimal && isValidIndependent && length(sel) < m) {
    return(list(verdict = "infeasible", witness = integer(0)))
  }
  list(verdict = "inconclusive", witness = integer(0))  # nocov
}

# ----- main -----------------------------------------------------------------

#' Exact Max-Min Diversity (MMDP) optimum on small instances
#'
#' Solves the Max-Min Diversity Problem (discrete p-dispersion) to proven
#' optimality by iterated node-packing \insertCite{Sayyady2016;textual}{MaxMin}: the optimum is
#' the largest threshold `lambda`, over the achieved distinct pairwise
#' distances, for which the threshold graph `G(lambda)` (edges join pairs
#' closer than `lambda`) contains an independent set of size at least `m`.
#' A binary search over the sorted distances resolves that threshold; each
#' probe solves a maximum-independent-set integer program with the `highs`
#' MILP backend. The problem is NP-hard, so this is intended only as an
#' external ground-truth reference on small instances, not a scalable method.
#'
#' @param d A `dist` object or a square symmetric numeric distance matrix.
#' @param m Integer target subset size, `2 <= m <= nrow(d)`.
#' @param solver Solver to use. Currently only `"highs"` is implemented;
#'   `NULL` selects it. Other values raise an error.
#' @param timeBudgetS Wall-clock budget in seconds for the whole search
#'   (shared across all internal IP solves). If the budget expires before the
#'   optimum is proven, the largest threshold proven feasible so far is
#'   returned with `proven = FALSE`.
#' @param progress Logical; show a progress bar during the binary search.
#'   Default: `TRUE` in interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
#' @return A list with fields
#'   \describe{
#'     \item{indices}{Integer vector of length `m`, sorted ascending: the
#'       selected points.}
#'     \item{objective}{The achieved `T_k` -- the minimum pairwise distance
#'       within `indices`. When `proven` is `TRUE` this equals the threshold
#'       `lambda` and is the true optimum; otherwise it is a valid lower
#'       bound on the optimum.}
#'     \item{proven}{Logical: `TRUE` if the search certified optimality
#'       within the budget, `FALSE` if it returned an unproven incumbent.}
#'     \item{time_s}{Wall-clock seconds elapsed.}
#'     \item{solver}{Name of the MILP backend used.}
#'     \item{n, m}{Instance size and target subset size.}
#'   }
#' @references \insertAllCited{}
#' @export
ExactMaxMin <- function(d, m, solver = NULL, timeBudgetS = 60,
                        progress = getOption("MaxMin.progress", interactive())) {
  t0 <- proc.time()[[3L]]
  if (is.null(solver)) solver <- "highs"
  if (!identical(solver, "highs")) {
    stop("Unsupported `solver`: ", solver, ". Only \"highs\" is implemented.")
  }
  if (!requireNamespace("highs", quietly = TRUE)) { # nocov start
    stop("The `highs` package is required for ExactMaxMin(). ",
         "Install it with install.packages(\"highs\").")
  } # nocov end

  d <- .ExactAsMatrix(d)
  n <- nrow(d)
  m <- as.integer(m)
  if (is.na(m) || m < 2L || m > n) {
    stop("`m` must satisfy 2 <= m <= nrow(d)")
  }

  Elapsed <- function() proc.time()[[3L]] - t0

  # Candidate thresholds: the achieved distinct pairwise distances, ascending.
  # The optimum is necessarily one of these (it is a realised distance), so
  # searching this finite grid is exact. Index k tests lambda = cand[k]; for
  # k = 1 (the smallest distance) the graph has no edges, so feasibility
  # holds whenever m <= n -- our guaranteed lower bound, established without
  # any IP solve.
  cand <- sort(unique(d[upper.tri(d)]))
  nCand <- length(cand)

  if (progress && nCand > 1L) {
    .pb <- cli::cli_progress_bar("ExactMaxMin", total = nCand - 1L,
                                 .auto_close = FALSE)
  }

  # Helper to package a result for a proven-feasible candidate index.
  Recover <- function(witness, lambda, proven) {
    idx <- sort(witness[seq_len(m)])
    sub <- d[idx, idx]
    diag(sub) <- Inf
    obj <- min(sub)
    # Tripwire (proven branch only): at the optimum the achieved minimum
    # equals lambda exactly for ANY m points drawn from the witness. A
    # mismatch would signal a solver or construction bug, never normal
    # degeneracy. On an unproven incumbent `obj` may exceed `lambda`; that
    # is a legitimately better lower bound, not an error.
    if (proven && abs(obj - lambda) > 1e-9 * max(1, abs(lambda))) { # nocov start
      stop("Internal error: recovered min-distance ", obj,
           " != proven threshold ", lambda)
    } # nocov end
    list(
      indices   = idx,
      objective = obj,
      proven    = proven,
      time_s    = Elapsed(),
      solver    = solver,
      n         = n,
      m         = as.integer(m)
    )
  }

  # Binary search for the largest feasible candidate index. Feasibility is
  # monotone non-increasing in lambda (raising the threshold only adds
  # edges, never removes them, so the independence number cannot grow), so
  # the feasible indices form a prefix [1 .. best] of the candidate vector.
  # Invariant: cand[1] is always feasible (empty graph, m <= n), so best
  # starts at 1 and we never need to prove the bottom of the range.
  lo <- 2L
  hi <- nCand
  bestIdx <- 1L
  bestWitness <- seq_len(n)            # any m points are >= cand[1] apart
  inconclusive <- FALSE

  while (lo <= hi) {
    remaining <- timeBudgetS - Elapsed()
    if (remaining <= 0) { # nocov start
      inconclusive <- TRUE
      break
    } # nocov end
    mid <- (lo + hi) %/% 2L
    v <- .MaxISVerdict(d, n, cand[mid], m, remaining)
    if (identical(v$verdict, "feasible")) {
      bestIdx <- mid
      bestWitness <- v$witness
      lo <- mid + 1L
      if (progress && nCand > 1L) {
        cli::cli_progress_update(id = .pb, set = (lo - 2L) + (nCand - hi))
      }
    } else if (identical(v$verdict, "infeasible")) {
      hi <- mid - 1L
      if (progress && nCand > 1L) {
        cli::cli_progress_update(id = .pb, set = (lo - 2L) + (nCand - hi))
      }
    } else { # nocov start
      # Budget expired mid-solve: cannot place this candidate. Stop and
      # return the best threshold proven feasible so far.
      inconclusive <- TRUE
      break
    } # nocov end
  }

  if (progress && nCand > 1L) {
    cli::cli_progress_done(id = .pb)
  }

  # Optimality is proven iff the search closed the interval (lo > hi) without
  # an inconclusive break: every candidate above bestIdx was certified
  # infeasible and bestIdx itself certified feasible.
  proven <- !inconclusive
  Recover(bestWitness, cand[bestIdx], proven)
}
