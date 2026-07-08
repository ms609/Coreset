# dropadd.R
#
# DropAdd Tabu Search (DropAdd-TS) for the Max-Min Diversity Problem.
#
# Implementation of Porumbel, Hao & Glover (2011) "A simple and effective
# algorithm for the MaxMin diversity problem", Annals of Operations Research
# 186:275-293. Algorithms 1-4 of that paper.
#
# The MMDPo (eq. 2 of the paper) objective is
#     max  [ min_{x,y in X} d(x,y) + eps * sum_{x,y in X} d(x,y) ]
# with a small eps. We use eps = 1e-9. The first term is the canonical MMDP
# (MaxMin) objective; the sum acts purely as a tie-break.
#
# Used in the manuscript as the SOTA heuristic competitor for MMDP per Marti
# et al. 2022 EJOR.

# ----- helpers --------------------------------------------------------------

# Constructive phase (Algorithm 1).
# Returns list(S = integer(k), iter_add = integer(k)) where iter_add[i] is
# the iteration at which S[i] was added (= i).
.DropAddConstruct <- function(dmat, k) {
  n <- nrow(dmat)
  S <- integer(k)
  # Seed point: argmax over Z of sum_y d(x, y) (Porumbel's eq. before Alg. 1).
  rowSumsD <- rowSums(dmat)
  S[1L] <- which.max(rowSumsD)  # ties: which.max → smallest index (deterministic)
  inS <- logical(n)
  inS[S[1L]] <- TRUE

  # Streamlined records.
  # minDist(x):       min over y in current_S of d(x, y), excluding x itself
  # sumDist(x):       sum over y in current_S of d(x, y) (d(x,x)=0 so harmless)
  # minDistCount(x): | { y in current_S, y != x : d(x, y) == minDist(x) } |
  minDist       <- dmat[, S[1L]]            # length n; minDist(S[1])=0 (self)
  minDist[S[1L]] <- Inf                     # never want to compare a selected
                                             # point to itself in min calcs
  sumDist       <- dmat[, S[1L]]
  minDistCount <- rep(1L, n)
  minDistCount[S[1L]] <- 0L                # no "other selected" yet

  if (k >= 2L) {
    for (h in 2L:k) {
      # Among AddX = Z \ S, pick x maximising (minDist, sumDist) lex,
      # ties → smallest index.
      cand <- which(!inS)
      md <- minDist[cand]
      bestMin <- max(md)
      tied <- cand[md == bestMin]
      if (length(tied) > 1L) {
        sd <- sumDist[tied]
        bestSum <- max(sd)
        tied <- tied[sd == bestSum]
      }
      xNew <- tied[1L]   # smallest index among lex-ties
      S[h] <- xNew
      inS[xNew] <- TRUE

      # Update streamlined records for ADD (Algorithm 3).
      # For every x != xNew:
      dCol <- dmat[, xNew]
      others <- seq_len(n)[-xNew]
      dToNew <- dCol[others]
      # sumDist
      sumDist[others] <- sumDist[others] + dToNew
      # Cases
      mdO <- minDist[others]
      caseB <- dToNew < mdO
      caseA <- dToNew == mdO & !caseB
      # Apply
      if (any(caseB)) {
        idx <- others[caseB]
        minDist[idx] <- dToNew[caseB]
        minDistCount[idx] <- 1L
      }
      if (any(caseA)) {
        idx <- others[caseA]
        minDistCount[idx] <- minDistCount[idx] + 1L
      }
      # xNew itself: its minDist becomes min over (S - {xNew}) of d.
      # Since we'd already maintained minDist[xNew] = Inf, we need to set
      # it now to the min distance from xNew to the previously-selected
      # points (which is what its minDist would have been without the
      # Inf override).
      prevS <- S[seq_len(h - 1L)]
      dXnew <- dmat[xNew, prevS]
      mn <- min(dXnew)
      minDist[xNew] <- mn
      minDistCount[xNew] <- sum(dXnew == mn)
    }
  }

  list(S = S, inS = inS, minDist = minDist, sumDist = sumDist,
       minDistCount = minDistCount)
}

# Internal test scaffolding (not part of the public API). Runs the production
# C++ DropAdd with tracing enabled and returns the dropped/added index
# sequences alongside the result, so the FIFO + tabu invariants can be asserted
# without exposing a `.trace` argument on `DropAdd()`. Kept deliberately thin:
# it mirrors only the `d`-path coercion and the single C++ call.
.DropAddTrace <- function(d, k, maxIter = NULL, plateau = 5000L,
                          maxSeconds = Inf) {
  dmat <- .AsDistMatrix(d)
  cppMaxIter <- if (is.null(maxIter)) .Machine$integer.max else as.integer(maxIter)
  out <- DropAdd_cpp(dmat, as.integer(k), as.double(maxSeconds),
                       cppMaxIter, as.integer(plateau), TRUE)
  list(
    indices = sort(as.integer(out$indices)),
    score   = as.numeric(out$objective),
    iters   = as.integer(out$iters),
    drops   = as.integer(out$drops),
    adds    = as.integer(out$adds)
  )
}

# ----- main entry point -----------------------------------------------------

#' DropAdd Tabu Search for the Max-Min Diversity Problem
#'
#' `DropAdd()` selects a maximally-dispersed subset of `k` points using the
#' DropAdd tabu search algorithm, which comprises a greedy construction
#' followed by a first-in, first-out drop-add tabu search, with streamlined
#' neighbour-evaluation tricks
#' \insertCite{@algorithms 1--4 in @Porumbel2011}{MaxMin}.
#'
#' @param k Integer; subset size, \eqn{2 \le k \le N}.
#' @param d A \code{dist} object or square symmetric numeric matrix.
#' @param points A numeric \eqn{N \times \mathrm{dim}} coordinate matrix (or an
#'  object coercible to one via \code{as.matrix}).
#'  Must be complete (no \code{NA}).
#'  Ignored if `d` specified.
#'  Avoids creating an \eqn{N \times N} distance matrix, enabling use at
#'  \eqn{N \ge 46340}).
#' @param plateau Integer; stop after this many consecutive drop-add
#'  iterations do not improve the score.
#' @param maxSeconds Numeric: terminate search after this many seconds have
#' elapsed.
#' @param seed Optional integer: a 1-based start index that overrides the
#'  construction's default warm-start seed (max-row-sum on the `d` path,
#'  centroid-peripheral on the `points` path), mirroring [FarFirst()]'s integer
#'  `strategy`. `NULL` (default) keeps the method's own seed. Not supported when
#'  candidate thinning binds (pass `maxCandidates = 0L`).
#' @templateVar default `46340L`
#' @templateVar default_basis the dense-distance-matrix feasibility ceiling
#'   (`floor(sqrt(.Machine$integer.max))`) the `points` path already crosses;
#'   below it nothing changes, at or above it the candidates are thinned to this
#'   size on the matrix-free path (no \eqn{m \times m} matrix is built)
#' @template maxCandidates
#' @templateVar progress_shows status messages are shown
#' @template progress
#'
#' @return `DropAdd()` returns an integer vector of length \code{k} containing the 1-based selected
#'   indices **sorted ascending** (unlike [FarFirst()], which returns
#'   farthest-first order), with attributes:
#'   \describe{
#'     \item{score}{numeric(1), achieved MaxMin objective
#'       \eqn{\min_{i \ne j \in S} d_{ij}}.}
#'     \item{secondary}{numeric(1), achieved sum of pairwise distances over
#'       \eqn{S} (upper-triangle sum).}
#'     \item{time_s}{numeric(1), wall-clock seconds spent.}
#'     \item{iters}{integer(1), main-loop iterations executed (excluding the
#'       construction phase).}
#'   }
#'   The vector has class `"MaxMinSelection"` and prints as a one-line summary
#'   (see [print.MaxMinSelection()]); it is otherwise an ordinary integer vector.
#'
#' @references \insertAllCited{}
#'
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(200), ncol = 2)
#' DropAdd(5L, dist(pts))
#'
#' # Composable coreset: thin to 40 candidates with farthest-first, then run
#' # DropAdd on the coreset. Returned indices are original-space row indices.
#' suppressWarnings(DropAdd(5L, points = pts, maxCandidates = 40L))
#'
#' # Disable thinning (run on the full problem):
#' DropAdd(5L, points = pts, maxCandidates = 0L)
#' @export
DropAdd <- function(k, d = NULL, plateau = 5000L, maxSeconds = Inf,
                    points = NULL, maxCandidates = 46340L, seed = NULL) {
  progress <- getOption("MaxMin.progress", interactive())
  if (!is.null(points) && !is.null(d)) {
    stop("supply `d` or `points`, not both")
  }
  usePoints <- !is.null(points)
  if (usePoints) {
    points <- .AsPointsMatrix(points)
    n <- nrow(points)
  } else {
    dmat <- .AsDistMatrix(d)
    n <- nrow(dmat)
  }
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > n) {
    stop("`k` must be a single integer with 2 <= k <= n")
  }
  plateau <- as.integer(plateau)
  if (length(plateau) != 1L || is.na(plateau) ||
      plateau < 1L) {
    stop("`plateau` must be a single positive integer")
  }
  if (!is.numeric(maxSeconds) || length(maxSeconds) != 1L ||
      is.na(maxSeconds) || maxSeconds <= 0) {
    stop("`maxSeconds` must be a single positive numeric (or Inf)")
  }
  # `seed` overrides the construction's default warm-start (max-row-sum on the
  # matrix path, centroid-peripheral on the points path) with an explicit 1-based
  # start index, mirroring FarFirst()'s integer `strategy`. NULL keeps the default.
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.finite(seed) || seed < 1L || seed > n)) {
    stop("`seed` must be a single index in [1, n], or NULL for the default")
  }
  seed0 <- if (is.null(seed)) -1L else as.integer(seed) - 1L

  # Composable-coreset thinning: when `maxCandidates` binds, farthest-first
  # reduces the n candidates to an m-point coreset and DropAdd runs on that.
  # The recursive call passes `maxCandidates = 0L`, so the coreset is solved by
  # the real kernel exactly once; on the points source the subproblem stays on
  # the matrix-free path (no m x m matrix is built).
  mc <- .ResolveCap(maxCandidates, n, k)
  if (!is.na(mc)) {
    if (!is.null(seed)) {
      stop("`seed=` is not supported with candidate thinning; pass ",
           "`maxCandidates = 0L` to run on the full problem")
    }
    return(.FarFirstThin(
      k, mc,
      d      = if (usePoints) NULL else dmat,
      points = if (usePoints) points else NULL,
      RunOnSubset = function(d, points)
        DropAdd(k, d = d, points = points, plateau = plateau,
                maxSeconds = maxSeconds, maxCandidates = 0L),
      label = "DropAdd"))
  }

  t0 <- proc.time()[[3L]]
  eps <- 1e-9

  # --- Matrix-free (coordinate) path ----------------------------------------
  if (usePoints) {
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (n = {n}, k = {k}, budget = {maxSeconds}s)",
        .auto_close = FALSE
      )
    }
    out <- DropAdd_points_cpp(points, k, as.double(maxSeconds),
                                .Machine$integer.max, plateau, FALSE, seed0)
    timeS <- proc.time()[[3L]] - t0
    if (progress) {
      itersMsg <- as.integer(out$iters)
      tkMsg    <- as.numeric(out$objective)
      cli::cli_process_done(
        msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
      )
    }
    return(.AsMaxMinSelection(structure(
      sort(as.integer(out$indices)),
      score     = as.numeric(out$objective),
      secondary = as.numeric(out$secondary),
      time_s    = timeS,
      iters     = as.integer(out$iters)
    ), "DropAdd"))
  }

  # --- C++ drop-add tabu search (the sole compute path) -------------------
  if (progress) {
    cli::cli_process_start(
      "DropAdd tabu search (n = {n}, k = {k}, budget = {maxSeconds}s)",
      .auto_close = FALSE
    )
  }
  out <- DropAdd_cpp(dmat, k, as.double(maxSeconds),
                       .Machine$integer.max, plateau, FALSE, seed0)
  timeS <- proc.time()[[3L]] - t0
  if (progress) {
    itersMsg <- as.integer(out$iters)
    tkMsg    <- as.numeric(out$objective)
    cli::cli_process_done(
      msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
    )
  }
  .AsMaxMinSelection(structure(
    sort(as.integer(out$indices)),
    score     = as.numeric(out$objective),
    secondary = as.numeric(out$secondary),
    time_s    = timeS,
    iters     = as.integer(out$iters)
  ), "DropAdd")
}
