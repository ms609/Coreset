# grasp.R
#
# GRASP with Path Relinking for the Max-Min Diversity Problem
# (Resende, Marti, Gallego & Duarte 2010, Computers & OR 37:498-508).
# Implements the static variant of Fig. 4.
#
# Public entry point: Grasp(). Internal helpers carry the .Grasp* prefix.
#
# This is a dense-matrix-only refinement metaheuristic: every phase
# (construction, local search, path relinking) operates on a materialised
# n x n distance matrix, so it does not offer a coordinate or column-oracle
# path and is intended for instances small enough to hold that matrix. For
# the matrix-free scaling regime use DropAdd(points=) / FarFirst(points=) or
# FarFirst() with a distance-column oracle.


# -- objective and selection helpers --------------------------------------

#' Minimum pairwise distance over a selection.
#' @keywords internal
.GraspObjective <- function(d, sel) {
  if (length(sel) < 2L) return(NA_real_)  # nocov
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  min(sub)
}

#' For each i in sel, return its nearest-other-selected distance.
#' @keywords internal
.GraspNearestInSel <- function(d, sel) {
  k <- length(sel)
  if (k < 2L) return(rep(Inf, k))  # nocov
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  apply(sub, 1L, min)
}

#' Count pairs at the minimum distance (used by extended-improvement LS).
#' @keywords internal
.GraspMinPairCount <- function(d, sel, dstar) {
  if (length(sel) < 2L) return(0L)  # nocov
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  # Each min-pair appears twice (i,j) and (j,i) in the symmetric matrix.
  sum(sub <= dstar) %/% 2L
}


# -- Phase 1: GRASP construction (S3) ------------------------------------

#' One randomised greedy construction.
#'
#' @param d Square distance matrix.
#' @param k Target subset size.
#' @param alpha RCL threshold parameter (alpha=1 -> greedy, alpha=0 -> random).
#' @return `.GraspConstruct()` returns an integer vector of length k.
#' @keywords internal
.GraspConstruct <- function(d, k, alpha) {
  n <- nrow(d)
  sel <- integer(k)
  sel[1L] <- sample.int(n, 1L)
  # nearest-selected distance for each candidate; Inf if no candidate available
  g <- d[, sel[1L]]
  g[sel[1L]] <- -Inf  # mark selected
  # `2L:k` counts *down* when k == 1L (-> c(2L, 1L)); guard so the helper is
  # safe when called directly with k == 1L, matching the C++ `h < k` loop.
  if (k >= 2L) {
    for (h in 2L:k) {
      # candidates are points with g > -Inf (i.e., not yet selected)
      candIdx <- which(g > -Inf)
      gv <- g[candIdx]
      gmax <- max(gv)
      gmin <- min(gv)
      thresh <- gmin + alpha * (gmax - gmin)
      rcl <- candIdx[gv >= thresh]
      # FP rounding at the documented alpha = 1 (and deterministically for
      # alpha > 1) can push thresh just past gmax, emptying the RCL. Fall back
      # to the unique greedy-best (the argmax-g candidate, first on ties)
      # WITHOUT an RNG draw, mirroring the size-1 branch and the C++ kernel so
      # R and C++ stay bit-identical.
      if (length(rcl) == 0L) {
        rcl <- candIdx[which.max(gv)]
      }
      if (length(rcl) == 1L) {
        pick <- rcl
      } else {
        pick <- rcl[sample.int(length(rcl), 1L)]
      }
      sel[h] <- pick
      # update g: candidate's new "min-to-sel" is min(old g, d to new pick)
      g <- pmin(g, d[, pick])
      g[pick] <- -Inf
    }
  }
  sort(sel)
}


# -- Phase 2: local search (S3.1, extended improvement) ------------------

#' Fast local search with the extended-improvement criterion.
#'
#' Iterates 1-swap moves on critical elements (those participating in a
#' min-distance edge). A swap is accepted if it strictly increases d*, or if
#' it preserves d* while reducing the count of pairs at d*.
#'
#' @param d Square distance matrix.
#' @param sel Integer vector of size k.
#' @return `.GraspLocalSearch()` returns an improved integer vector of size k.
#' @keywords internal
.GraspLocalSearch <- function(d, sel) {
  n <- nrow(d)
  k <- length(sel)
  if (k < 2L) return(sel)  # nocov
  repeat {
    selSorted <- sort(sel)
    di <- .GraspNearestInSel(d, selSorted)
    dstar <- min(di)
    pairCount <- .GraspMinPairCount(d, selSorted, dstar)
    critical <- which(di <= dstar)  # indices into selSorted
    inSel <- logical(n); inSel[selSorted] <- TRUE
    outIdx <- which(!inSel)
    bestSwap <- NULL
    bestDstar <- dstar
    bestPaircount <- pairCount
    for (ci in critical) {
      drop <- selSorted[ci]
      remaining <- selSorted[-ci]
      for (s in outIdx) {
        cand <- c(remaining, s)
        # quick: compute d* directly
        sub <- d[cand, cand, drop = FALSE]
        diag(sub) <- Inf
        newDstar <- min(sub)
        if (newDstar > bestDstar) {
          bestSwap <- list(drop = drop, add = s)
          bestDstar <- newDstar
          bestPaircount <- sum(sub <= newDstar) %/% 2L
        } else if (newDstar == bestDstar) {
          newPaircount <- sum(sub <= newDstar) %/% 2L
          if (newPaircount < bestPaircount) {
            bestSwap <- list(drop = drop, add = s)
            bestPaircount <- newPaircount
          }
        }
      }
    }
    if (is.null(bestSwap)) break
    sel <- c(selSorted[selSorted != bestSwap$drop], bestSwap$add)
  }
  sort(sel)
}


# -- Phase 3: path relinking (S4.1, greedy) ------------------------------

#' Greedy path relinking from x toward y.
#'
#' @param d Square distance matrix.
#' @param x,y Integer selections of equal length.
#' @return `.GraspPathRelink()` returns `list(best = best selection on path,
#'   intermediates = number of intermediate states visited including endpoints)`.
#' @keywords internal
.GraspPathRelink <- function(d, x, y) {
  x <- sort(x); y <- sort(y)
  if (identical(x, y)) {
    return(list(best = x, objective = .GraspObjective(d, x), intermediates = 0L))
  }
  toDrop <- setdiff(x, y)  # |toDrop| = r
  toAdd  <- setdiff(y, x)
  r <- length(toDrop)
  pk <- x
  bestSel <- x
  bestZ <- .GraspObjective(d, x)
  zY <- .GraspObjective(d, y)
  if (zY > bestZ) { bestSel <- y; bestZ <- zY }
  intermediates <- 0L
  for (k in seq_len(r)) {
    # pk_minus_y = toDrop intersect pk; pk_plus_y = toAdd \ pk
    dropCands <- intersect(pk, toDrop)
    addCands  <- setdiff(toAdd, pk)
    bestPair <- NULL
    bestPairZ <- -Inf
    # tie-breaker: lexicographic on (drop_index, add_index) -> iterate sorted
    for (i in sort(dropCands)) {
      for (j in sort(addCands)) {
        cand <- c(pk[pk != i], j)
        sub <- d[cand, cand, drop = FALSE]
        diag(sub) <- Inf
        zc <- min(sub)
        if (zc > bestPairZ) {
          bestPairZ <- zc
          bestPair <- c(i, j)
        }
      }
    }
    pk <- sort(c(pk[pk != bestPair[1L]], bestPair[2L]))
    intermediates <- intermediates + 1L
    if (bestPairZ > bestZ) {
      bestSel <- pk
      bestZ <- bestPairZ
    }
  }
  list(best = bestSel, objective = bestZ, intermediates = intermediates)
}


# -- Elite-set helpers ---------------------------------------------------

#' Hamming distance in selection space: m - |intersection|.
#' @keywords internal
.GraspHammingToES <- function(sel, ES) {
  m <- length(sel)
  vapply(ES, function(e) m - length(intersect(sel, e)), integer(1L))
}

#' Try to insert sel into the elite set ES.
#'
#' @return `.GraspTryInsert()` returns the updated ES (list of selections, sorted best-to-worst by z).
#' @keywords internal
.GraspTryInsert <- function(d, ES, esZ, sel, selZ, dth) {
  z1 <- esZ[1L]
  zb <- esZ[length(esZ)]
  hamm <- .GraspHammingToES(sel, ES)
  dmin <- min(hamm)
  accept <- FALSE
  if (selZ > z1) accept <- TRUE
  else if (selZ > zb && dmin >= dth) accept <- TRUE
  # Don't admit a duplicate
  if (dmin == 0L) accept <- FALSE
  if (!accept) return(list(ES = ES, esZ = esZ, changed = FALSE))
  # Eviction is restricted to members WORSE than the candidate -- Resende et al.
  # (2010) Fig. 4 line 8, "closest solution to x' in ES with z(x') > z(x^k)", and
  # §4.1, "we remove the closest solution to x' in ES among those worse than it
  # in value".
  #
  # The pool is never empty: acceptance already requires selZ > z1 (every member
  # is then worse) or selZ > zb (the worst member qualifies).
  worse <- which(esZ < selZ)
  hammWorse <- hamm[worse]
  cand <- worse[hammWorse == min(hammWorse)]
  # Ties: lowest z, then lowest index (which.min gives the first minimum).
  closest <- cand[which.min(esZ[cand])]
  ES <- ES[-closest]
  esZ <- esZ[-closest]
  # insert sorted (descending z)
  pos <- sum(esZ >= selZ) + 1L
  if (pos > length(ES)) {
    ES <- c(ES, list(sel))
    esZ <- c(esZ, selZ)
  } else {
    ES <- c(ES[seq_len(pos - 1L)], list(sel), ES[pos:length(ES)])
    esZ <- c(esZ[seq_len(pos - 1L)], selZ, esZ[pos:length(esZ)])
  }
  list(ES = ES, esZ = esZ, changed = TRUE)
}


# -- Outer loop: GRASP with PR static variant (Fig. 4) -------------------

#' GRASP with Path Relinking for the Max-Min Diversity Problem
#'
#' `Grasp()` solves the Max-Min Diversity Problem (discrete _p_-dispersion) with
#' the static variant of the GRASP / path-relinking metaheuristic
#' \insertCite{@@@@Resende2010, fig. 4}{MaxMin}{}. This is the most expensive
#' heuristic in this package, and attains correspondingly high-quality
#' selections.
#'
#' The GRASP with path-relinking algorithm conducts a randomised-greedy
#' construction with extended-improvement local search builds; it identifies an
#' elite set, then conducts a single pass of path relinking over all elite
#' pairs \insertCite{Resende2010}{MaxMin}.
#'
#' The refinement loop stops after `plateau` consecutive GRASP iterations
#' fail to improve the best elite objective, or once `maxSeconds` have
#' elapsed.
#'
#' This method will fail if the complete \eqn{N \times N} distance matrix is
#' too large to fit into memory.
#'
#' @param k Integer subset size, `2 <= k <= nrow(d)`.
#' @param d Either a `dist` object or a square symmetric numeric matrix.
#' @param plateau Integer; stop after this many consecutive GRASP
#'   iterations without an improvement to the best elite objective. The
#'   primary, deterministic stopping criterion.
#' @param eliteSize Size of the elite set |ES|.
#' @param alpha Construction greediness in `[0, 1]`. Each step draws the next
#'   point at random from a shortlist of the strongest candidates -- those
#'   whose gain lies within a fraction `alpha` of the best-to-worst spread.
#'   `alpha = 1` is pure greedy (best only); `alpha = 0` is uniform random
#'   among candidates.
#' @param maxSeconds Numeric specifying wall-clock ceiling, in seconds.
#' @templateVar default `2000L`
#' @templateVar default_basis conservative because `Grasp()` is matrix-only, so
#'   the coreset subproblem is a dense \eqn{m \times m} matrix
#'   (\eqn{2000 \times 2000 \approx} 32 MB)
#' @template maxCandidates
#' @return `Grasp()` returns an integer vector of length `k` specifying the
#' indices of the selected points, with attributes:
#'   \describe{
#'     \item{score}{Achieved MaxMin objective \eqn{T_k}.}
#'     \item{time_s}{Wall-clock seconds spent.}
#'     \item{iters}{Number of GRASP refinement iterations executed.}
#'     \item{pr_calls}{Number of path-relinking pair-applications run.}
#'   }
#'   The vector has class `"MaxMinSelection"` and prints as a one-line summary
#'   (see [print.MaxMinSelection()]).
#' @templateVar progress_shows a bar tracks how close the search is to its `plateau` stopping criterion, snapping back each time a better solution is found
#' @template progress
#' @references \insertAllCited{}
#'
#' @seealso [DropAdd()] for scalable refinement;
#'   [ExactMaxMin()] for the proven optimum on small instances.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' Grasp(5L, dist(pts), plateau = 20L, eliteSize = 4L)
#'
#' # Composable coreset: thin to 20 candidates with farthest-first, then run
#' # GRASP on the coreset. Returned indices are original-space row indices.
#' suppressWarnings(Grasp(5L, dist(pts), plateau = 20L, maxCandidates = 20L))
#' @export
Grasp <- function(k, d, plateau = 100L, eliteSize = 10L, alpha = 0.8,
                  maxSeconds = Inf, maxCandidates = 2000L) {
  progress <- getOption("MaxMin.progress", interactive())
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  k <- as.integer(k)
  stopifnot(k >= 2L, k <= n)
  eliteSize <- as.integer(eliteSize)
  stopifnot(eliteSize >= 1L)
  plateau <- as.integer(plateau)
  stopifnot(plateau >= 1L)
  if (!is.numeric(maxSeconds) || length(maxSeconds) != 1L ||
      is.na(maxSeconds) || maxSeconds <= 0) {
    stop("`maxSeconds` must be a single positive numeric (or Inf)")
  }
  # `alpha` outside [0, 1] is meaningless for the RCL threshold and, for
  # alpha > 1, deterministically empties the RCL (see .GraspConstruct); reject
  # it before it reaches the kernel.
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number in [0, 1]")
  }

  # Composable-coreset thinning: when `maxCandidates` binds, farthest-first
  # reduces the n candidates to an m-point coreset (a cheap d[core, core] subset
  # of the matrix Grasp already holds) and Grasp runs on that. The recursive
  # call passes `maxCandidates = 0L`, so the coreset is solved exactly once.
  mc <- .ResolveCap(maxCandidates, n, k)
  if (!is.na(mc)) {
    return(.FarFirstThin(
      k, mc, d = d, points = NULL,
      RunOnSubset = function(d, points)
        Grasp(k, d = d, plateau = plateau, eliteSize = eliteSize,
              alpha = alpha, maxSeconds = maxSeconds, maxCandidates = 0L),
      label = "Grasp"))
  }

  # The kernel reports its stall counter (consecutive non-improving GRASP
  # iterations, 0..plateau) through `cb`; we render it as a bar that fills as the
  # search stagnates toward the `plateau` stopping criterion and snaps back to 0
  # each time a better solution is found. cb is read-only, so passing it cannot
  # change the result. Errors in the bar must never abort a long search.
  cb <- NULL
  if (progress) {
    pb <- cli::cli_progress_bar(
      format = paste0(
        "{cli::pb_spin} GRASP-PR (n = {n}, k = {k}) | ",
        "stall {cli::pb_current}/{cli::pb_total} {cli::pb_bar} {cli::pb_percent}"
      ),
      total = plateau, clear = FALSE, .auto_close = FALSE
    )
    cb <- function(noImprove) {
      tryCatch(
        cli::cli_progress_update(set = noImprove, id = pb, force = TRUE),
        error = function(e) NULL
      )
    }
  }

  out <- Grasp_cpp(d, k, plateau, .Machine$integer.max, eliteSize,
                     as.double(alpha), as.double(maxSeconds), cb)

  if (progress) {
    cli::cli_progress_done(id = pb)
    itersMsg <- as.integer(out$iters)
    tkMsg    <- as.numeric(out$objective)
    timeMsg  <- as.numeric(out$time_s)
    cli::cli_alert_success(
      "Grasp: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeMsg, 1)}s"
    )
  }

  # Return:
  .AsMaxMinSelection(structure(
    sort(as.integer(out$indices)),
    score    = as.numeric(out$objective),
    time_s   = as.numeric(out$time_s),
    iters    = as.integer(out$iters),
    pr_calls = as.integer(out$pr_calls)
  ), "Grasp")
}

# Pure-R reference implementation of Grasp, used as the parity oracle for
# the compiled kernel (see tests/testthat/test-grasp.R). It mirrors
# Grasp_cpp() step for step, including the construction RNG draws, so that
# from a common `set.seed()` the two agree bit for bit. Not exported; callers
# use Grasp().
#
# Termination is the deterministic stagnation rule: stop after
# `plateau` consecutive iterations that do not raise the best elite
# objective esZ[1] (which is monotone non-decreasing under .GraspTryInsert()).
# `maxIter` is an optional hard cap; `maxSeconds` an optional ceiling
# (Inf = off) that leaves the result reproducible.
#' @keywords internal
.Grasp_R <- function(k, d, plateau, maxIter = .Machine$integer.max,
                       eliteSize = 10L, alpha = 0.8, maxSeconds = Inf) {
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  k <- as.integer(k)
  eliteSize <- as.integer(eliteSize)
  dth <- 5L
  t0 <- proc.time()[[3L]]

  # Phase A: build initial elite set.
  ES <- vector("list", eliteSize)
  esZ <- numeric(eliteSize)
  for (b in seq_len(eliteSize)) {
    x  <- .GraspConstruct(d, k, alpha)
    xp <- .GraspLocalSearch(d, x)
    ES[[b]] <- xp
    esZ[b] <- .GraspObjective(d, xp)
  }
  ord <- order(esZ, decreasing = TRUE)
  ES <- ES[ord]; esZ <- esZ[ord]

  bestSel <- ES[[1L]]
  bestZ   <- esZ[1L]
  iters   <- 0L
  prCalls <- 0L
  on.exit(setTimeLimit(), add = TRUE)
  # The wall-clock ceiling is the only stopping criterion that survives when
  # both `plateau` and `maxIter` are disabled, so it must actually fire.
  # setTimeLimit() alone cannot be trusted for this: a non-positive `elapsed`
  # *removes* the limit, so once Phase A has already overrun a tiny budget,
  # `maxSeconds - elapsed <= 0` would silently disable it and Phase B would
  # spin forever. Guard with an explicit predicate checked each iteration, and
  # only arm setTimeLimit() -- which interrupts a single overlong iteration --
  # while time genuinely remains.
  budgetSpent <- function() {
    is.finite(maxSeconds) && (proc.time()[[3L]] - t0) >= maxSeconds
  }
  if (is.finite(maxSeconds) && !budgetSpent()) {
    setTimeLimit(elapsed = maxSeconds - (proc.time()[[3L]] - t0),
                 transient = TRUE)
  }

  # Phase B: GRASP iterations until `plateau` consecutive non-improving
  # iterations (the deterministic criterion), an optional iteration cap, or an
  # optional wall-clock ceiling.
  tryCatch({
    noImprove <- 0L
    repeat {
      if (budgetSpent()) break
      if (noImprove >= plateau) break
      if (iters >= maxIter) break
      x  <- .GraspConstruct(d, k, alpha)
      xp <- .GraspLocalSearch(d, x)
      zp <- .GraspObjective(d, xp)
      res <- .GraspTryInsert(d, ES, esZ, xp, zp, dth)
      ES <- res$ES; esZ <- res$esZ
      iters <- iters + 1L
      if (esZ[1L] > bestZ) {
        bestZ   <- esZ[1L]
        bestSel <- ES[[1L]]
        noImprove <- 0L
      } else {
        noImprove <- noImprove + 1L
      }
    }

    # Phase C: path relinking over all elite pairs (deterministic; no RNG).
    bestSel <- ES[[1L]]
    bestZ   <- esZ[1L]
    k <- length(ES)
    if (k >= 2L) {
      for (i in seq_len(k - 1L)) {
        for (j in (i + 1L):k) {
          pr1 <- .GraspPathRelink(d, ES[[i]], ES[[j]])
          pr2 <- .GraspPathRelink(d, ES[[j]], ES[[i]])
          prCalls <- prCalls + 2L
          ySel <- if (pr1$objective >= pr2$objective) pr1$best else pr2$best
          yp <- .GraspLocalSearch(d, ySel)
          zp <- .GraspObjective(d, yp)
          if (zp > bestZ) {
            bestZ <- zp
            bestSel <- yp
          }
        }
      }
    }
  }, error = function(e) { # nocov start
    setTimeLimit()
    if (!grepl("time limit", conditionMessage(e), ignore.case = TRUE)) stop(e)
  }) # nocov end

  # Return:
  .AsMaxMinSelection(structure(
    sort(as.integer(bestSel)),
    score    = bestZ,
    time_s   = proc.time()[[3L]] - t0,
    iters    = iters,
    pr_calls = prCalls
  ), "Grasp")
}
