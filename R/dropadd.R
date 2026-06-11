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
# Returns list(S = integer(m), iter_add = integer(m)) where iter_add[k] is
# the iteration at which S[k] was added (= k).
.DropAddConstruct <- function(dmat, m) {
  n <- nrow(dmat)
  S <- integer(m)
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

  if (m >= 2L) {
    for (h in 2L:m) {
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

# nocov start
# MMDPo objective from streamlined records: for the current set X,
#   maxmin = min over x in X of minDist(x)
.MMDPoObj <- function(dmat, S) {
  if (length(S) < 2L) return(c(maxmin = NA_real_, sumpair = 0))
  sub <- dmat[S, S, drop = FALSE]
  diag(sub) <- Inf
  mm <- min(sub)
  sp <- sum(sub[is.finite(sub)]) / 2  # half because symmetric
  # use upper tri instead for clarity:
  diag(sub) <- 0
  sp <- sum(sub[upper.tri(sub)])
  c(maxmin = mm, sumpair = sp)
}

# Recompute minDist(x) and minDistCount(x) from scratch given current S
# (Alg. 4 recompute branch).
.RecomputeMinFor <- function(x, dmat, S) {
  # Excludes x itself if it happens to be in S.
  others <- if (x %in% S) setdiff(S, x) else S
  if (!length(others)) return(c(Inf, 0L))
  vals <- dmat[x, others]
  mn <- min(vals)
  list(mn = mn, cnt = sum(vals == mn))
}
# nocov end

# ----- main entry point -----------------------------------------------------

#' DropAdd Tabu Search for the Max-Min Diversity Problem
#'
#' Implements the DropAdd-TS algorithm of \insertCite{Porumbel2011;textual}{MaxMin} for
#' selecting a maximally-dispersed subset of \code{m} points from a distance
#' matrix. The procedure consists of a deterministic greedy construction
#' (Algorithm 1) followed by a FIFO drop-add tabu search (Algorithm 2) with
#' the streamlined neighbour-evaluation tricks of Algorithms 3 and 4.
#'
#' The MMDPo objective optimised is
#' \deqn{\min_{x,y \in X} d(x,y) + \epsilon \sum_{x,y \in X} d(x,y),}
#' with \eqn{\epsilon = 10^{-9}}, so the pairwise sum acts as a tie-break.
#'
#' Tabu mechanics. The algorithm maintains an integer \code{iter} stamp for
#' every point, set to the iteration at which the point last entered or
#' left the selected set \eqn{X}. At each main-loop iteration the point with
#' the smallest \code{iter} value (the oldest member, FIFO) is dropped, and
#' the point in \eqn{Z \setminus X} maximising lexicographically
#' \eqn{(\mathrm{minDist}, \mathrm{sumDist})} is added. The FIFO
#' invariant guarantees that across any window of \eqn{m} iterations every
#' initially-selected point is dropped exactly once before any re-eviction.
#'
#' Time budget behaviour. The time budget (\code{timeBudgetS}) is checked at
#' most once every 256 iterations (matrix-free path) or 1024 iterations
#' (matrix path). On large instances where each iteration is slow, the actual
#' elapsed time may exceed the specified budget by up to one iteration's worth
#' of computation.
#'
#' @param d A \code{dist} object or square symmetric numeric matrix.
#'   Mutually exclusive with \code{points}; supply exactly one.
#' @param points A numeric \eqn{n \times \mathrm{dim}} coordinate matrix (or an
#'   object coercible to one via \code{as.matrix}). Must be complete (no
#'   \code{NA}). Mutually exclusive with \code{d}; supply exactly one. When
#'   supplied, the algorithm never materialises the dense \eqn{n \times n}
#'   distance matrix, giving \eqn{O(n)} working memory and enabling use at
#'   \eqn{n} far exceeding the matrix path's ceiling (R's
#'   \code{as.matrix.dist} overflows at \eqn{n = 46340}).
#' @param m Integer; subset size, \eqn{2 \le m \le n}.
#' @param plateau Integer; stop after this many consecutive drop-add
#'   iterations that do not improve the best objective. The primary,
#'   deterministic stopping criterion. The search is RNG-free (ties broken by
#'   smallest index), so for a given instance the result is reproducible and
#'   machine-independent. Default 5000.
#' @param maxIter Optional integer hard cap on iterations (excluding
#'   construction). \code{NULL} (default) leaves \code{plateau} in sole
#'   control.
#' @param timeBudgetS Optional wall-clock ceiling in seconds, checked at
#'   iteration boundaries. Default \code{Inf} (no ceiling, fully reproducible).
#'   A finite value caps runtime but makes the result machine-dependent.
#' @param progress Logical; show a start/done status line. Default: `TRUE` in
#'   interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`). No effect when
#'   `.verify = TRUE` (testing path).
#' @param .verify Logical (testing only); if `TRUE`, routes to the R reference
#'   loop and brute-force asserts the streamlined records at every iteration.
#'   Default `FALSE` (the C++ fast path). Only applies to the \code{d} path.
#' @param .trace Optional environment (testing only); if supplied, the dropped
#'   and added index sequences are written into it as `drops` and `adds`.
#'
#' @return An integer vector of length \code{m} containing the 1-based selected
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
#'   (see [print.MaxMinSelection]); it is otherwise an ordinary integer vector.
#'
#' @references \insertAllCited{}
#'
#' @export
DropAdd <- function(d = NULL, m, plateau = 5000L, maxIter = NULL,
                      timeBudgetS = Inf,
                      progress = getOption("MaxMin.progress", interactive()),
                      points = NULL,
                      .verify = FALSE, .trace = NULL) {
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
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L || m > n) {
    stop("`m` must be a single integer with 2 <= m <= n")
  }
  plateau <- as.integer(plateau)
  if (length(plateau) != 1L || is.na(plateau) ||
      plateau < 1L) {
    stop("`plateau` must be a single positive integer")
  }
  if (!is.null(maxIter)) {
    maxIter <- as.integer(maxIter)
    if (length(maxIter) != 1L || is.na(maxIter) || maxIter < 0L) {
      stop("`maxIter` must be NULL or a single non-negative integer")
    }
  }
  if (!is.numeric(timeBudgetS) || length(timeBudgetS) != 1L ||
      is.na(timeBudgetS) || timeBudgetS <= 0) {
    stop("`timeBudgetS` must be a single positive numeric (or Inf)")
  }

  t0 <- proc.time()[[3L]]
  eps <- 1e-9

  # --- Matrix-free (coordinate) path ----------------------------------------
  if (usePoints) {
    cppMaxIter <- if (is.null(maxIter)) .Machine$integer.max else maxIter
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (n = {n}, m = {m}, budget = {timeBudgetS}s)",
        .auto_close = FALSE
      )
    }
    out <- DropAdd_points_cpp(points, m, as.double(timeBudgetS),
                                cppMaxIter, plateau, FALSE)
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

  # --- C++ fast path. .verify routes to R for the brute-force assertion. --
  if (!.verify) {
    cppMaxIter <- if (is.null(maxIter)) .Machine$integer.max else maxIter
    wantTrace <- !is.null(.trace)
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (n = {n}, m = {m}, budget = {timeBudgetS}s)",
        .auto_close = FALSE
      )
    }
    out <- DropAdd_cpp(dmat, m, as.double(timeBudgetS),
                         cppMaxIter, plateau, wantTrace)
    if (wantTrace) {
      .trace$drops <- out$drops
      .trace$adds  <- out$adds
    }
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

  # --- Construction (Algorithm 1) -----------------------------------------
  cons <- .DropAddConstruct(dmat, m)
  S            <- cons$S
  inS          <- cons$inS
  minDist      <- cons$minDist
  sumDist      <- cons$sumDist
  minDistCount <- cons$minDistCount

  # Track best-known MMDPo solution. Compute the initial objective from the
  # streamlined records (not via .MMDPoObj's upper-triangle sum) so that the
  # R fallback path and the C++ port are bit-identical from the first iter:
  # `.MMDPoObj` sums `dmat[S, S]` upper-tri in row-major order, the
  # streamlined `sum(sumDist[S]) / 2` accumulates m partial sums — both
  # compute the same quantity but differ in the last ULP.
  bestS <- S
  bestMaxmin  <- min(minDist[S])
  bestSumpair <- sum(sumDist[S]) / 2
  bestScore   <- bestMaxmin + eps * bestSumpair

  # --- Drop-Add tabu search (Algorithm 2) ---------------------------------
  # FIFO via circular buffer: S is a length-m vector treated as a queue.
  # `head` indexes the oldest member. Each iteration drops S[head], writes
  # the new x* into S[head], and advances head modulo m. The Porumbel
  # iter_stamp / which.min(iter_stamp[S]) pair is redundant under this
  # invariant — head IS the FIFO head.
  head        <- 1L
  itersDone   <- 0L
  noImprove   <- 0L
  effectiveMax <- if (is.null(maxIter)) .Machine$integer.max else maxIter
  if (m >= n) effectiveMax <- 0L  # all points selected: no drop-add move exists
  if (!is.null(.trace)) {
    .trace$drops <- integer(0)
    .trace$adds  <- integer(0)
  }

  # Internal brute-force assertion of streamlined records.
  .VerifyRecords <- function(S, minDist, minDistCount, sumDist) {
    for (xx in seq_len(n)) {
      others <- if (xx %in% S) setdiff(S, xx) else S
      if (!length(others)) next  # nocov
      vals <- dmat[xx, others]
      tm <- min(vals); tc <- sum(vals == tm); ts <- sum(dmat[xx, S])
      stopifnot(abs(minDist[xx] - tm) < 1e-9,
                minDistCount[xx] == tc,
                abs(sumDist[xx] - ts) < 1e-6)
    }
    invisible(TRUE)
  }
  if (.verify) .VerifyRecords(S, minDist, minDistCount, sumDist)

  on.exit(setTimeLimit(), add = TRUE)
  # The budget is the only stopping rule that survives when both `effectiveMax`
  # and `plateau` are disabled, so it must actually fire. setTimeLimit() alone
  # is unsafe: a non-positive `elapsed` *removes* the limit, so a budget already
  # spent during construction would silently disable it and the loop would spin
  # forever. Arm setTimeLimit() (to interrupt a single overlong iteration) only
  # while time remains, and re-check the budget at the foot of the loop so it
  # always halts -- after at least one iteration -- regardless.
  budgetSpent <- function() {
    is.finite(timeBudgetS) && (proc.time()[[3L]] - t0) >= timeBudgetS
  }
  if (is.finite(timeBudgetS) && !budgetSpent()) {
    setTimeLimit(elapsed = timeBudgetS - (proc.time()[[3L]] - t0),
                 transient = TRUE)
  }
  tryCatch(repeat {
    if (itersDone >= effectiveMax) break
    if (noImprove >= plateau) break

    # 1. DROP: head IS the FIFO head under the circular-buffer invariant.
    xHash <- S[head]

    # Apply drop to streamlined records (Algorithm 4). Full-vector updates
    # operate on all n entries; dmat[xHash, xHash] = 0 so sumDist[xHash]
    # is unchanged, and xHash's minDist/minDistCount get clobbered by
    # the trailing self-recompute regardless.
    inS[xHash] <- FALSE
    dDropFull <- dmat[, xHash]
    sumDist <- sumDist - dDropFull
    affectedMask <- dDropFull <= minDist
    affectedMask[xHash] <- FALSE         # xHash gets its own recompute below
    if (any(affectedMask)) {
      minDistCount[affectedMask] <- minDistCount[affectedMask] - 1L
    }
    needRecompute <- which(affectedMask & minDistCount == 0L)
    sAfterDrop <- S[-head]
    if (length(needRecompute)) {
      # Vectorised recompute. `sub` is k × (m-1); rows for needRecompute
      # entries still in S need their self entry masked to Inf.
      sub <- dmat[needRecompute, sAfterDrop, drop = FALSE]
      inSRows <- inS[needRecompute]
      if (any(inSRows)) {
        sel <- which(inSRows)
        cols <- match(needRecompute[sel], sAfterDrop)
        sub[cbind(sel, cols)] <- Inf
      }
      mns  <- as.numeric(do.call(pmin.int, asplit(sub, 2L)))
      cnts <- rowSums(sub == mns)
      finiteMns <- is.finite(mns)
      if (all(finiteMns)) {
        minDist[needRecompute]       <- mns
        minDistCount[needRecompute] <- cnts
      } else {
        bad <- !finiteMns
        minDist[needRecompute[bad]]       <- Inf
        minDistCount[needRecompute[bad]] <- 0L
        if (any(finiteMns)) {
          good <- needRecompute[finiteMns]
          minDist[good]       <- mns[finiteMns]
          minDistCount[good] <- cnts[finiteMns]
        }
      }
    }
    # xHash's own minDist record: distance to nearest surviving peer.
    if (length(sAfterDrop)) {
      rowVals <- dmat[xHash, sAfterDrop]
      mn <- min(rowVals)
      minDist[xHash] <- mn
      minDistCount[xHash] <- sum(rowVals == mn)
    } else {
      minDist[xHash] <- Inf         # nocov
      minDistCount[xHash] <- 0L   # nocov
    }

    # 2. ADD: argmax over Add X(k) = Z - X(k) of (minDist, sumDist), ties →
    # smallest idx. xHash is excluded for this iteration (Porumbel et al. 2011,
    # p.281): the just-dropped point cannot be re-added immediately — the tabu
    # rule that prevents looping. It is eligible again once head advances.
    cand <- which(!inS)
    cand <- cand[cand != xHash]
    md <- minDist[cand]
    bestMin <- max(md)
    tied <- cand[md == bestMin]
    if (length(tied) > 1L) {
      sd <- sumDist[tied]
      bestSum <- max(sd)
      tied <- tied[sd == bestSum]
    }
    xNew <- tied[1L]

    # Update streamlined records for ADD (Algorithm 3). Same full-vector
    # discipline; xNew's own record is overwritten by the trailing block.
    inS[xNew] <- TRUE
    dColFull <- dmat[, xNew]
    sumDist <- sumDist + dColFull
    caseB <- dColFull <  minDist
    caseA <- dColFull == minDist        # < and == are disjoint; no !caseB mask needed
    if (any(caseB)) {
      minDist[caseB]       <- dColFull[caseB]
      minDistCount[caseB] <- 1L
    }
    if (any(caseA)) {
      minDistCount[caseA] <- minDistCount[caseA] + 1L
    }
    # xNew's own minDist: distance to peers in S \ {xHash}.
    rowVals2 <- dmat[xNew, sAfterDrop]
    mn2 <- min(rowVals2)
    minDist[xNew] <- mn2
    minDistCount[xNew] <- sum(rowVals2 == mn2)

    # Write xNew into the head slot and advance.
    S[head] <- xNew
    head <- if (head == m) 1L else head + 1L

    # 3. Test for improvement of best-known MMDPo solution.
    curMaxmin  <- min(minDist[S])
    curSumpair <- sum(sumDist[S]) / 2     # double-counted
    curScore   <- curMaxmin + eps * curSumpair
    if (curScore > bestScore) {
      bestS       <- S
      bestMaxmin  <- curMaxmin
      bestSumpair <- curSumpair
      bestScore   <- curScore
      noImprove   <- 0L
    } else {
      noImprove   <- noImprove + 1L
    }

    itersDone <- itersDone + 1L

    if (!is.null(.trace)) {
      .trace$drops <- c(.trace$drops, xHash)
      .trace$adds  <- c(.trace$adds, xNew)
    }
    if (.verify) .VerifyRecords(S, minDist, minDistCount, sumDist)
    if (budgetSpent()) break
  }, error = function(e) { # nocov start
    setTimeLimit()
    if (!grepl("time limit", conditionMessage(e), ignore.case = TRUE)) stop(e)
  }) # nocov end

  timeS <- proc.time()[[3L]] - t0
  # Return:
  .AsMaxMinSelection(structure(
    sort(as.integer(bestS)),
    score     = as.numeric(bestMaxmin),
    secondary = as.numeric(bestSumpair),
    time_s    = timeS,
    iters     = as.integer(itersDone)
  ), "DropAdd")
}
