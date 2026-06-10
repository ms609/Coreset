# grasp.R
#
# GRASP with Path Relinking for the Max-Min Diversity Problem
# (Resende, Marti, Gallego & Duarte 2010, Computers & OR 37:498-508).
# Implements the static variant of Fig. 4.
#
# Public entry point: GraspPR(). Internal helpers carry the .Gpr* prefix.
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
.GprObjective <- function(d, sel) {
  if (length(sel) < 2L) return(NA_real_)  # nocov
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  min(sub)
}

#' For each i in sel, return its nearest-other-selected distance.
#' @keywords internal
.GprNearestInSel <- function(d, sel) {
  m <- length(sel)
  if (m < 2L) return(rep(Inf, m))  # nocov
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  apply(sub, 1L, min)
}

#' Count pairs at the minimum distance (used by extended-improvement LS).
#' @keywords internal
.GprMinPairCount <- function(d, sel, dstar) {
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
#' @param m Target subset size.
#' @param alpha RCL threshold parameter (alpha=1 -> greedy, alpha=0 -> random).
#' @return Integer vector of length m.
#' @keywords internal
.GprConstruct <- function(d, m, alpha) {
  n <- nrow(d)
  sel <- integer(m)
  sel[1L] <- sample.int(n, 1L)
  # nearest-selected distance for each candidate; Inf if no candidate available
  g <- d[, sel[1L]]
  g[sel[1L]] <- -Inf  # mark selected
  for (h in 2L:m) {
    candMask <- is.finite(g) & g > -Inf
    # candidates are points with g > -Inf (i.e., not yet selected)
    candIdx <- which(g > -Inf)
    gv <- g[candIdx]
    gmax <- max(gv)
    gmin <- min(gv)
    thresh <- gmin + alpha * (gmax - gmin)
    rcl <- candIdx[gv >= thresh]
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
#' @param sel Integer vector of size m.
#' @return Improved integer vector of size m.
#' @keywords internal
.GprLocalSearch <- function(d, sel) {
  n <- nrow(d)
  m <- length(sel)
  if (m < 2L) return(sel)  # nocov
  repeat {
    selSorted <- sort(sel)
    di <- .GprNearestInSel(d, selSorted)
    dstar <- min(di)
    pairCount <- .GprMinPairCount(d, selSorted, dstar)
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
#' @return list(best = best selection on path, intermediates = number of
#'   intermediate states visited including endpoints).
#' @keywords internal
.GprPathRelink <- function(d, x, y) {
  x <- sort(x); y <- sort(y)
  if (identical(x, y)) {
    return(list(best = x, objective = .GprObjective(d, x), intermediates = 0L))
  }
  toDrop <- setdiff(x, y)  # |toDrop| = r
  toAdd  <- setdiff(y, x)
  r <- length(toDrop)
  pk <- x
  bestSel <- x
  bestZ <- .GprObjective(d, x)
  zY <- .GprObjective(d, y)
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
.GprHammingToES <- function(sel, ES) {
  m <- length(sel)
  vapply(ES, function(e) m - length(intersect(sel, e)), integer(1L))
}

#' Try to insert sel into the elite set ES.
#'
#' @return Updated ES (list of selections, sorted best-to-worst by z).
#' @keywords internal
.GprTryInsert <- function(d, ES, esZ, sel, selZ, dth) {
  z1 <- esZ[1L]
  zb <- esZ[length(esZ)]
  hamm <- .GprHammingToES(sel, ES)
  dmin <- min(hamm)
  accept <- FALSE
  if (selZ > z1) accept <- TRUE
  else if (selZ > zb && dmin >= dth) accept <- TRUE
  # Don't admit a duplicate
  if (dmin == 0L) accept <- FALSE
  if (!accept) return(list(ES = ES, esZ = esZ, changed = FALSE))
  # remove the closest member (smallest Hamming distance); break ties on lowest z
  closest <- which(hamm == dmin)
  if (length(closest) > 1L) {
    closest <- closest[which.min(esZ[closest])]
  } else {
    closest <- closest[1L]
  }
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
#' Solves the Max-Min Diversity Problem (discrete p-dispersion) with the
#' GRASP / path-relinking metaheuristic of \insertCite{Resende2010;textual}{MaxMin},
#' static variant (their Fig. 4): a randomised-greedy construction
#' with extended-improvement local search builds and maintains an elite set,
#' followed by a single pass of path relinking over all elite pairs. On the
#' application benchmark this attains the highest \eqn{T_k} of the methods in
#' this package, at correspondingly higher cost.
#'
#' **Deterministic termination.** The refinement loop stops after
#' `plateau` consecutive GRASP iterations that fail to improve the best
#' elite objective (rather than after a wall-clock budget). Given a fixed
#' `seed`, the entire run — construction RNG, iteration count, and result — is
#' therefore reproducible and machine-independent; `set.seed(seed)` controls
#' the same random stream the compiled kernel consumes via R's RNG. An
#' optional `timeBudgetS` ceiling is available as a safety cap, but using a
#' finite value reintroduces machine-dependence and is off by default.
#'
#' This is a **dense-matrix-only** method: it materialises and repeatedly
#' subsets the full \eqn{n \times n} distance matrix, so it is suited to
#' instances small enough to hold that matrix. It offers no coordinate or
#' column-oracle path. For the matrix-free regime where the dense matrix is
#' infeasible, use [DropAdd()] (coordinate path via \code{points =}) or [FarFirst()] (coordinate or
#' distance-column oracle path), whose
#' \eqn{T_k} lands within roughly a percent on the benchmark while scaling to
#' far larger instances.
#'
#' @param d Either a `dist` object or a square symmetric numeric matrix.
#' @param m Integer subset size, `2 <= m <= nrow(d)`.
#' @param plateau Integer; stop after this many consecutive GRASP
#'   iterations without an improvement to the best elite objective. The
#'   primary, deterministic stopping criterion. Default 100.
#' @param maxIter Optional integer hard cap on GRASP refinement iterations
#'   (excluding the elite-set construction). `NULL` (default) leaves
#'   `plateau` in sole control.
#' @param eliteSize Size of the elite set |ES|. Default 10.
#' @param alpha RCL threshold; `alpha = 1` is pure greedy, `alpha = 0`
#'   uniform random. Default 0.8.
#' @param timeBudgetS Optional wall-clock ceiling in seconds. Default `Inf`
#'   (no ceiling, fully reproducible). A finite value caps runtime but makes
#'   the result machine-dependent.
#' @param seed Optional integer; if supplied, `set.seed(seed)` is called at
#'   entry. `GraspPR` is genuinely stochastic (randomised construction and RCL
#'   sampling), so the seed governs the trajectory and the returned selection.
#' @return An integer vector of length `m` (1-based, sorted ascending)
#'   with attributes:
#'   \describe{
#'     \item{score}{Achieved MaxMin objective \eqn{T_k}.}
#'     \item{time_s}{Wall-clock seconds spent.}
#'     \item{iters}{Number of GRASP refinement iterations executed.}
#'     \item{pr_calls}{Number of path-relinking pair-applications run.}
#'   }
#' @references \insertAllCited{}
#'
#' @seealso [DropAdd()] for scalable refinement;
#'   [ExactMaxMin()] for the proven optimum on small instances.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' res <- GraspPR(dist(pts), m = 5L, plateau = 20L, eliteSize = 4L,
#'                seed = 1L)
#' res
#' @export
GraspPR <- function(d, m, plateau = 100L, maxIter = NULL,
                    eliteSize = 10L, alpha = 0.8, timeBudgetS = Inf,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  m <- as.integer(m)
  stopifnot(m >= 2L, m <= n)
  eliteSize <- as.integer(eliteSize)
  stopifnot(eliteSize >= 1L)
  plateau <- as.integer(plateau)
  stopifnot(plateau >= 1L)
  if (is.null(maxIter)) {
    maxIter <- .Machine$integer.max
  } else {
    maxIter <- as.integer(maxIter)
    stopifnot(maxIter >= 0L)
  }
  if (!is.numeric(timeBudgetS) || length(timeBudgetS) != 1L ||
      is.na(timeBudgetS) || timeBudgetS <= 0) {
    stop("`timeBudgetS` must be a single positive numeric (or Inf)")
  }

  out <- GraspPR_cpp(d, m, plateau, maxIter, eliteSize,
                     as.double(alpha), as.double(timeBudgetS))
  structure(
    sort(as.integer(out$indices)),
    score    = as.numeric(out$objective),
    time_s   = as.numeric(out$time_s),
    iters    = as.integer(out$iters),
    pr_calls = as.integer(out$pr_calls)
  )
}

# Pure-R reference implementation of GraspPR, used as the parity oracle for
# the compiled kernel (see tests/testthat/test-gpr.R). It mirrors
# GraspPR_cpp() step for step, including the construction RNG draws, so that
# from a common `set.seed()` the two agree bit for bit. Not exported; callers
# use GraspPR().
#
# Termination is the deterministic stagnation rule: stop after
# `plateau` consecutive iterations that do not raise the best elite
# objective esZ[1] (which is monotone non-decreasing under .GprTryInsert()).
# `maxIter` is an optional hard cap; `timeBudgetS` an optional ceiling
# (Inf = off) that leaves the result reproducible.
#' @keywords internal
.GraspPR_R <- function(d, m, plateau, maxIter = .Machine$integer.max,
                       eliteSize = 10L, alpha = 0.8, timeBudgetS = Inf) {
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  m <- as.integer(m)
  eliteSize <- as.integer(eliteSize)
  dth <- 5L
  t0 <- proc.time()[[3L]]

  # Phase A: build initial elite set.
  ES <- vector("list", eliteSize)
  esZ <- numeric(eliteSize)
  for (b in seq_len(eliteSize)) {
    x  <- .GprConstruct(d, m, alpha)
    xp <- .GprLocalSearch(d, x)
    ES[[b]] <- xp
    esZ[b] <- .GprObjective(d, xp)
  }
  ord <- order(esZ, decreasing = TRUE)
  ES <- ES[ord]; esZ <- esZ[ord]

  bestSel <- ES[[1L]]
  bestZ   <- esZ[1L]
  iters   <- 0L
  prCalls <- 0L
  on.exit(setTimeLimit(), add = TRUE)
  if (is.finite(timeBudgetS)) {
    setTimeLimit(elapsed = max(0, timeBudgetS - (proc.time()[[3L]] - t0)),
                 transient = TRUE)
  }

  # Phase B: GRASP iterations until `plateau` consecutive non-improving
  # iterations (the deterministic criterion), an optional iteration cap, or an
  # optional wall-clock ceiling.
  tryCatch({
    noImprove <- 0L
    repeat {
      if (noImprove >= plateau) break
      if (iters >= maxIter) break
      x  <- .GprConstruct(d, m, alpha)
      xp <- .GprLocalSearch(d, x)
      zp <- .GprObjective(d, xp)
      res <- .GprTryInsert(d, ES, esZ, xp, zp, dth)
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
          pr1 <- .GprPathRelink(d, ES[[i]], ES[[j]])
          pr2 <- .GprPathRelink(d, ES[[j]], ES[[i]])
          prCalls <- prCalls + 2L
          ySel <- if (pr1$objective >= pr2$objective) pr1$best else pr2$best
          yp <- .GprLocalSearch(d, ySel)
          zp <- .GprObjective(d, yp)
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

  structure(
    sort(as.integer(bestSel)),
    score    = bestZ,
    time_s   = proc.time()[[3L]] - t0,
    iters    = iters,
    pr_calls = prCalls
  )
}
