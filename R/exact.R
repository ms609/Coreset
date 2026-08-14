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
# with min-distance >= lambda.
#
# Two structural choices make this practical well beyond toy sizes:
#
#   * Each probe is decided as a k-clique search on the complement graph H
#     (pairs >= lambda apart), by the branch-and-bound in exact_reduce.cpp:
#     peel to the (k-1)-core, split into components, and search each under a
#     greedy-colouring bound. The equivalent packing IP spends its time
#     closing an LP relaxation whose bound is far above k, and needs a
#     constraint matrix of one row per G-edge.
#
#   * Heuristic warm start + galloping search. Rather than bisecting the whole
#     sorted distance vector (~log2(n^2/2) probes, half of them splitting hairs
#     around the optimum), we seed a provably-achievable lower bound from the
#     package's own heuristics (Grasp restarts + DropAdd), then gallop upward
#     from there to the first infeasible threshold and bisect the small
#     bracket. When a heuristic attains the optimum -- common at the small k of
#     interest -- a single infeasibility proof certifies it. The seed only sets
#     the starting lower bound: the search still proves the true optimum
#     regardless of seed quality (a loose seed merely costs extra probes, never
#     a wrong answer).
#
# Used in the manuscript ONLY as an external ground-truth reference on small
# instances -- it is NP-hard and not a scalable method.

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
  # The warm start calls Grasp() and DropAdd() inside a tryCatch, so their own
  # symmetry guard would surface only as an empty heuristic pool and a silent
  # fall back to the trivial bound; the check belongs here. A non-finite entry
  # fails a symmetry test on its own terms (NA != NA), so it is left to the
  # guard that names it.
  if (AllFinite_cpp(d, .NThreads()) && !IsExactlySymmetric_cpp(d)) {
    stop("`d` must be symmetric: d[i, j] and d[j, i] must be equal. ",
         "Use `(d + t(d)) / 2` if rounding has made them differ.")
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
  # `maxCandidates = 0L` pins both to no thinning: the warm starts must see the
  # full matrix the user handed to the exact solver, regardless of the solvers'
  # own default coreset caps.
  grasps <- local({
    op <- options(MaxMin.progress = FALSE); on.exit(options(op))
    lapply(seq_len(nStart), function(s)
      tryCatch(Grasp(k, d, plateau = 50L, maxCandidates = 0L),
               error = function(e) NULL))
  })
  raw <- c(if (is.null(warmStart)) list() else list(warmStart),
           grasps,
           list(tryCatch(local({
             op <- options(MaxMin.progress = FALSE); on.exit(options(op))
             DropAdd(d = d, k = k, plateau = 512L, maxCandidates = 0L)
           }), error = function(e) NULL)))
  pool <- Filter(Negate(is.null), lapply(raw, valid))
  if (!length(pool)) {
    return(NULL)
  }
  vals <- vapply(pool, scv, numeric(1))
  j <- which.max(vals)
  list(value = vals[[j]], witness = pool[[j]])
}

# Decide one maximum-independent-set feasibility probe on the threshold graph
# G(lambda) against the target size k.
#
# `hi`, `hj` are the endpoints of the pairs at least lambda apart -- the edges
# of the complement graph H, in which an independent k-set of G(lambda) is a
# k-clique, and which `ThresholdDecide_cpp` searches exhaustively.
#
# A returned witness is validated against `d` independently of what produced
# it: the selected set is checked to contain no G(lambda) edge and its size is
# compared to k.
#
# Returns list(verdict, witness) where verdict is one of:
#   "feasible"     -- a validated independent set of size >= k was found
#                     (witness = its vertex indices; min-distance >= lambda),
#   "infeasible"   -- the search was exhaustive and no k-clique exists
#                     (witness = integer(0)),
#   "inconclusive" -- the budget expired before either could be established
#                     (witness = integer(0)).
.MaxISVerdict <- function(d, n, hi, hj, lambda, k, timeLimit) {
  if (!is.finite(timeLimit) || timeLimit <= 0) { # nocov start
    return(list(verdict = "inconclusive", witness = integer(0)))
  } # nocov end

  if (length(hi) == as.double(n) * (n - 1) / 2) {  # n^2 overflows an integer
    # H holds every pair, so G(lambda) is edgeless: every vertex is
    # independent and alpha = n. Feasible whenever k <= n (the caller's
    # guard).
    return(list(verdict = "feasible", witness = seq_len(n)))
  }

  res <- ThresholdDecide_cpp(hi, hj, n, k, timeLimit)
  witness <- res[["witness"]]
  if (identical(res[["status"]], "feasible")) {
    sub <- d[witness, witness, drop = FALSE]
    if (length(witness) != k || any(sub[upper.tri(sub)] < lambda)) { # nocov start
      stop("Internal error: clique search returned a set violating threshold ",
           lambda)
    } # nocov end
    return(list(verdict = "feasible", witness = witness))
  }
  list(verdict = res[["status"]], witness = integer(0))
}

# ----- main -----------------------------------------------------------------

#' Exact Max-Min Diversity Problem solution
#'
#' `ExactMaxMin()` finds the optimal solution to the Max-Min Diversity Problem
#' (discrete _p_-dispersion) by iterated node-packing
#' \insertCite{Sayyady2016}{MaxMin} (which may be slow or intractable on large
#' sets).
#'
#' The search is warm-started from a heuristic lower bound (the best of several
#' [Grasp()] restarts and a [DropAdd()] pass), then gallops upward from that
#' bound to the first infeasible threshold and bisects the resulting bracket.
#' When a heuristic already attains the optimum, a single infeasibility proof
#' certifies it.
#' Each feasibility probe is first reduced to its \eqn{(k-1)}-core and greedily
#' coloured, then searched exhaustively for a witness under a colouring bound.
#' The search runs on one core: its branches parallelise, but measurably only
#' for infeasibility proofs, and threads would make the reported subset
#' thread-dependent.
#' The indices returned may vary between releases where several subsets attain
#' the optimum; the `score` does not.
#'
#' @param k Integer: target subset size, between 2 and `nrow(d)`.
#' @param d `dist` object or a square symmetric numeric distance matrix.
#' @param maxSeconds Numeric: search terminates after this many seconds have
#' elapsed, returning largest threshold proven feasible.
#' @param warmStart Integer vector giving indices of a candidate subset to add
#'  to the heuristic warm-start pool, e.g. a selection computed by another
#'  solver.
#' @templateVar progress_shows a progress indicator is shown
#' @template progress
#' @return `ExactMaxMin()` returns an integer vector of length `k` (sorted
#'   ascending) with class `"MaxMinSelection"`, carrying attributes:
#'   \describe{
#'     \item{score}{The minimum pairwise distance within the selection.
#'       When `proven` is `TRUE` this is the optimum; otherwise
#'       a lower bound.}
#'     \item{proven}{Logical: `TRUE` if the search certified optimality within
#'       the budget, `FALSE` if it returned an unproven incumbent.}
#'     \item{time_s}{Wall-clock seconds elapsed.}
#'     \item{N, k}{Instance size and target subset size.}
#'   }
#'   Prints as a terse summary via [print.MaxMinSelection()].
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(18), ncol = 2)
#' ExactMaxMin(3L, dist(pts))
#' @export
ExactMaxMin <- function(k, d, maxSeconds = 60, warmStart = NULL) {
  progress <- getOption("MaxMin.progress", interactive())
  t0 <- proc.time()[[3L]]
  d <- .ExactAsMatrix(d)
  n <- nrow(d)
  k <- as.integer(k)
  if (is.na(k) || k < 2L || k > n) {
    stop("`k` must satisfy 2 <= k <= nrow(d)")
  }

  Elapsed <- function() proc.time()[[3L]] - t0

  # Warm start: a provably-achievable lower bound + witness. Every candidate
  # <= ws$value is feasible (the witness attains it), so the optimum is at
  # least ws$value. Without a heuristic, fall back to the trivial bound (the
  # smallest distance, whose threshold graph is edgeless).
  ws <- .ExactWarmStart(d, n, k, warmStart)
  lowest <- if (is.null(ws)) -Inf else ws$value # nocov

  # Candidate thresholds: the achieved distinct distances from `lowest` up,
  # ascending. The optimum is necessarily one of these (it is a realised
  # distance) and is at least `lowest`, so searching this finite grid is exact.
  # Radix sort then adjacent dedupe -- exactly the ascending distinct values.
  s <- sort.int(TriangleAtLeast_cpp(d, lowest), method = "radix")
  cand <- s[c(TRUE, s[-1L] != s[-length(s)])]
  nCand <- length(cand)

  if (progress) {
    .pb <- cli::cli_progress_bar("ExactMaxMin", total = NA, .auto_close = FALSE)
  }
  tick <- function() if (progress) cli::cli_progress_update(id = .pb) # nocov

  # Feasibility oracle at a candidate index: does G(cand[idx]) admit an
  # independent set of size >= k?
  feas <- function(idx, remaining) {
    lambda <- cand[idx]
    h <- EdgesAtLeast_cpp(d, lambda)
    .MaxISVerdict(d, n, h[["hi"]], h[["hj"]], lambda, k, remaining)
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
        N      = n,
        k      = as.integer(k)
      ),
      producer = "ExactMaxMin"
    )
  }

  # cand[1] is the warm start's own realised value, so the search starts there.
  i0 <- 1L
  bestIdx <- 1L
  if (is.null(ws)) { # nocov start
    bestWitness <- seq_len(n)
  } else { # nocov end
    bestWitness <- ws$witness
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
