# exact.R
#
# Exact solver for the Max-Min Diversity Problem (MMDP = discrete
# p-dispersion): given an n x n symmetric distance matrix and a target size
# k, find the k-subset S maximising T_k = min_{i != j in S} d(i, j).
#
# Method: iterated node-packing (Sayyady & Fathi 2016, EJOR 253(1):216-225).
# The MaxMin optimum is the largest threshold lambda over the achieved
# distinct pairwise distances such that the threshold graph
#   G(lambda) = (V, { (i, j) : d(i, j) < lambda })
# admits an independent set of size >= k. An independent set in G(lambda) is
# a set of points all pairwise >= lambda apart -- i.e. a feasible k-subset
# with min-distance >= lambda. At each tested lambda we solve a maximum-
# independent-set feasibility IP
#   maximize sum_i x_i  s.t.  x_i + x_j <= 1 for every edge (i, j);  x in {0,1}
# declaring lambda feasible iff the optimum (the independence number) >= k.
#
# Two structural choices make this practical well beyond toy sizes:
#
#   * Sparse packing matrix. The constraint matrix carries exactly two
#     non-zeros per edge; it is built as a sparse `dgCMatrix` rather than a
#     dense nEdge x n matrix. The dense form is O(nEdge * n) memory -- several
#     gigabytes per solve once n is a few hundred -- and is the reason a naive
#     node-packing solver stalls at a couple of hundred points.
#
#   * Heuristic warm start + galloping search. Rather than bisecting the whole
#     sorted distance vector (~log2(n^2/2) solves, half of them splitting hairs
#     around the optimum), we seed a provably-achievable lower bound from the
#     package's own heuristics (Grasp restarts + DropAdd), then gallop upward
#     from there to the first infeasible threshold and bisect the small
#     bracket. When a heuristic attains the optimum -- common at the small k of
#     interest -- a single infeasibility solve certifies it. The seed only sets
#     the starting lower bound: the search still proves the true optimum
#     regardless of seed quality (a loose seed merely costs extra solves, never
#     a wrong answer).
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

# A provably-achievable k-subset, as a lower bound on the optimum. Grasp
# (several RNG restarts) attains the package's best heuristic T_k on small-to-
# medium instances; DropAdd adds a deterministic anchor. The best-achieving
# subset wins. `warmStart`, if a valid k-subset, joins the pool. Returns
# list(value, witness), or NULL if no heuristic produced a usable subset.
.ExactWarmStart <- function(d, n, k, warmStart, nStart = 8L) {
  scv   <- function(idx) { s <- d[idx, idx]; diag(s) <- Inf; min(s) }
  # Coerce a raw selection to a valid sorted k-subset, or drop it (NULL).
  valid <- function(idx) {
    idx <- tryCatch(sort(unique(as.integer(idx))), error = function(e) integer(0))
    if (length(idx) == k && idx[1L] >= 1L && idx[k] <= n) idx else NULL
  }
  # Heuristic pool: optional caller warmStart, several Grasp restarts (drawing
  # on the session RNG -- this is where their diversity comes from), one
  # deterministic DropAdd. Like Grasp() itself, this advances the session RNG.
  grasps <- local({
    op <- options(MaxMin.progress = FALSE); on.exit(options(op))
    lapply(seq_len(nStart), function(s)
      tryCatch(Grasp(k, d, plateau = 50L), error = function(e) NULL))
  })
  raw <- c(if (is.null(warmStart)) list() else list(warmStart),
           grasps,
           list(tryCatch(local({
             op <- options(MaxMin.progress = FALSE); on.exit(options(op))
             DropAdd(d = d, k = k, plateau = 512L)
           }), error = function(e) NULL)))
  pool <- Filter(Negate(is.null), lapply(raw, valid))
  if (!length(pool)) {
    return(NULL)
  }
  vals <- vapply(pool, scv, numeric(1))
  j <- which.max(vals)
  list(value = vals[[j]], witness = pool[[j]])
}

# Solve one maximum-independent-set feasibility probe on the threshold graph
# G(lambda) and classify the result against the target size k.
#
# `ei`, `ej` are the endpoint index vectors of the edges (pairs with
# d(i, j) < lambda), supplied by the caller from a one-off upper-triangle
# precompute. The packing constraint x_i + x_j <= 1 is assembled as a SPARSE
# matrix (two non-zeros per edge), then sum x is maximised over binary
# variables. The returned witness is validated independently of the solver
# status: x is rounded, the selected set is checked to contain no G(lambda)
# edge, and its size is compared to k. This makes the classification robust to
# solver status-code or time-limit quirks.
#
# Returns list(verdict, witness) where verdict is one of:
#   "feasible"     -- a validated independent set of size >= k was found
#                     (witness = its vertex indices; min-distance >= lambda),
#   "infeasible"   -- the IP was solved to proven optimality and the maximum
#                     independent set has size < k (witness = integer(0)),
#   "inconclusive" -- the budget expired before either could be established
#                     (witness = integer(0)).
.MaxISVerdict <- function(d, n, ei, ej, lambda, k, timeLimit) {
  if (!is.finite(timeLimit) || timeLimit <= 0) { # nocov start
    return(list(verdict = "inconclusive", witness = integer(0)))
  } # nocov end
  nEdge <- length(ei)

  if (nEdge == 0L) {
    # Empty graph: every vertex is independent, so alpha = n. Feasible
    # whenever k <= n (guaranteed by the caller's guard). No solve needed.
    return(list(verdict = "feasible", witness = seq_len(n)))
  }

  # One row per edge: x_i + x_j <= 1. Sparse: exactly two non-zeros per row.
  A <- Matrix::sparseMatrix(
    i = rep.int(seq_len(nEdge), 2L),
    j = c(ei, ej),
    x = 1, dims = c(nEdge, n)
  )

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

  if (isValidIndependent && length(sel) >= k) {
    return(list(verdict = "feasible", witness = sel))
  }

  # Infeasibility is provable ONLY when the IP reached optimality and the
  # certified maximum independent set is still smaller than k. A time-limit
  # hit with a too-small (or empty) incumbent proves nothing.
  optimal <- identical(res$status_message, "Optimal")
  if (optimal && isValidIndependent && length(sel) < k) {
    return(list(verdict = "infeasible", witness = integer(0)))
  }
  list(verdict = "inconclusive", witness = integer(0))  # nocov
}

# ----- main -----------------------------------------------------------------

#' Exact Max-Min Diversity Problem solution
#'
#' `ExactMaxMin()` finds the optimal solution to the Max-Min Diversity Problem
#' (discrete _p_-dispersion) by iterated node-packing
#' \insertCite{Sayyady2016;textual}{MaxMin}.
#' As this problem is NP-hard, it is feasible only for small sets.
#'
#' The search is warm-started from a heuristic lower bound (the best of several
#' [Grasp()] restarts and a [DropAdd()] pass), then gallops upward from that
#' bound to the first infeasible threshold and bisects the resulting bracket.
#' When a heuristic already attains the optimum, a single infeasibility solve
#' certifies it.
#'
#' The proven `objective` is exact and does not depend on the RNG. Only the
#' returned `indices` can vary when several subsets attain the optimum: the warm
#' start draws on the session RNG via [Grasp()] and, like `Grasp()`, advances it.
#' Call [set.seed()] before `ExactMaxMin()` for a reproducible selection.
#'
#' @param k Integer: target subset size, between 2 and `nrow(d)`.
#' @param d `dist` object or a square symmetric numeric distance matrix.
#' @param maxSeconds Numeric: search terminates after this many seconds have
#' elapsed; returning largest threshold proven feasible so far.
#' @param warmStart Integer vector giving indices of a candidate subset to add
#'  to the heuristic warm-start pool, e.g. a selection computed by another
#'  solver.
#' @templateVar progress_shows a progress indicator is shown
#' @template progress
#' @return `ExactMaxMin()` returns an integer vector of length `k` (sorted
#'   ascending) with class `"MaxMinSelection"`, carrying attributes:
#'   \describe{
#'     \item{score}{Achieved `T_k` -- the minimum pairwise distance within the
#'       selection. When `proven` is `TRUE` this is the true optimum; otherwise
#'       a lower bound.}
#'     \item{proven}{Logical: `TRUE` if the search certified optimality within
#'       the budget, `FALSE` if it returned an unproven incumbent.}
#'     \item{time_s}{Wall-clock seconds elapsed.}
#'     \item{solver}{The name of the MILP backend: `"highs"`.}
#'     \item{N, k}{Instance size and target subset size.}
#'   }
#'   Prints as a one-line summary via [print.MaxMinSelection()].
#'   `inherits(result, "MaxMinSelection")` is `TRUE`.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(18), ncol = 2)
#' ExactMaxMin(3L, dist(pts))
#' @export
ExactMaxMin <- function(k, d, maxSeconds = 60, warmStart = NULL) {
  progress <- getOption("MaxMin.progress", interactive())
  t0 <- proc.time()[[3L]]
  if (!requireNamespace("highs", quietly = TRUE)) { # nocov start
    stop("The `highs` package is required for ExactMaxMin(). ",
         "Install it with install.packages(\"highs\").")
  } # nocov end

  d <- .ExactAsMatrix(d)
  n <- nrow(d)
  k <- as.integer(k)
  if (is.na(k) || k < 2L || k > n) {
    stop("`k` must satisfy 2 <= k <= nrow(d)")
  }

  Elapsed <- function() proc.time()[[3L]] - t0

  # Upper-triangle edge structure, precomputed ONCE: every probe is then a
  # threshold over `ud` (a vector compare) rather than a fresh n x n logical.
  ut <- which(upper.tri(d))
  rc <- arrayInd(ut, dim(d))
  ui <- rc[, 1L]; uj <- rc[, 2L]; ud <- d[ut]

  # Candidate thresholds: the achieved distinct pairwise distances, ascending.
  # The optimum is necessarily one of these (it is a realised distance), so
  # searching this finite grid is exact. cand[1] (the smallest distance) gives
  # an edgeless graph, feasible whenever k <= n -- the guaranteed lower bound,
  # established without any IP solve.
  cand <- sort(unique(ud))
  nCand <- length(cand)

  if (progress) {
    .pb <- cli::cli_progress_bar("ExactMaxMin", total = NA, .auto_close = FALSE)
  }
  tick <- function() if (progress) cli::cli_progress_update(id = .pb) # nocov

  # Feasibility oracle at a candidate index: does G(cand[idx]) admit an
  # independent set of size >= k?
  feas <- function(idx, remaining) {
    lambda <- cand[idx]
    e <- (ud < lambda)
    .MaxISVerdict(d, n, ui[e], uj[e], lambda, k, remaining)
  }

  # Helper to package a result for a proven-feasible candidate index.
  Recover <- function(witness, lambda, proven) {
    idx <- sort(witness[seq_len(k)])
    sub <- d[idx, idx]
    diag(sub) <- Inf
    obj <- min(sub)
    # Tripwire (proven branch only): at the optimum the achieved minimum
    # equals lambda exactly -- the largest feasible threshold has the next
    # candidate proven infeasible, so the witness cannot beat it. A mismatch
    # would signal a construction bug, never normal degeneracy. On
    # an unproven incumbent `obj` may exceed `lambda`; that is a legitimately
    # better lower bound, not an error.
    if (proven && abs(obj - lambda) > 1e-9 * max(1, abs(lambda))) { # nocov start
      stop("Internal error: recovered min-distance ", obj,
           " != proven threshold ", lambda)
    } # nocov end
    .AsMaxMinSelection(
      structure(
        idx,
        score  = obj,
        proven = proven,
        time_s = Elapsed(),
        solver = "highs",
        N      = n,
        k      = as.integer(k)
      ),
      producer = "ExactMaxMin"
    )
  }

  # Warm start: a provably-achievable lower bound + witness. Every candidate
  # <= ws$value is feasible (the witness attains it), so the optimum index is
  # at least i0. Without a heuristic, fall back to the trivial bound (cand[1]).
  ws <- .ExactWarmStart(d, n, k, warmStart)
  if (is.null(ws)) { # nocov start
    i0 <- 1L; bestIdx <- 1L; bestWitness <- seq_len(n)
  } else { # nocov end
    i0 <- findInterval(ws$value, cand)        # cand[i0] == ws$value (realised)
    bestIdx <- i0; bestWitness <- ws$witness
  }
  inconclusive <- FALSE

  # Gallop up from i0 to the first infeasible threshold. Feasibility is
  # monotone (raising lambda only adds edges, never grows the independence
  # number), so the feasible indices form a prefix [1 .. best]; doubling the
  # step finds the boundary in O(log gap) when the warm start is near-optimal.
  loF <- i0; hiX <- NA_integer_; step <- 1L; probe <- i0 + 1L
  while (probe <= nCand) {
    rem <- maxSeconds - Elapsed()
    if (rem <= 0) { inconclusive <- TRUE; break } # nocov
    v <- feas(probe, rem); tick()
    if (identical(v$verdict, "feasible")) {
      loF <- probe; bestIdx <- probe; bestWitness <- v$witness
      step <- step * 2L; probe <- probe + step
    } else if (identical(v$verdict, "infeasible")) {
      hiX <- probe; break
    } else { # nocov start
      inconclusive <- TRUE; break
    } # nocov end
  }
  if (is.na(hiX)) hiX <- nCand + 1L           # nothing above proven infeasible

  # Bisect the bracket (loF, hiX): the largest feasible index in between.
  if (!inconclusive) {
    lo <- loF + 1L
    hi <- min(hiX - 1L, nCand)
    while (lo <= hi) {
      rem <- maxSeconds - Elapsed()
      if (rem <= 0) { inconclusive <- TRUE; break } # nocov
      mid <- (lo + hi) %/% 2L
      v <- feas(mid, rem); tick()
      if (identical(v$verdict, "feasible")) {
        bestIdx <- mid; bestWitness <- v$witness; lo <- mid + 1L
      } else if (identical(v$verdict, "infeasible")) {
        hi <- mid - 1L
      } else { # nocov start
        inconclusive <- TRUE; break
      } # nocov end
    }
  }

  if (progress) cli::cli_progress_done(id = .pb) # nocov

  # Optimality is proven iff the search closed without an inconclusive break:
  # bestIdx certified feasible (heuristic witness or IP witness) and the next
  # candidate certified infeasible (or bestIdx is the largest distance).
  proven <- !inconclusive
  Recover(bestWitness, cand[bestIdx], proven)
}
