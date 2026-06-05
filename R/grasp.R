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
# the matrix-free scaling regime use DropAddTSPoints() / Gonzalez(points=) or
# Gonzalez() with a distance-column oracle.


# -- objective and selection helpers --------------------------------------

#' Minimum pairwise distance over a selection.
#' @keywords internal
.GprObjective <- function(d, sel) {
  if (length(sel) < 2L) return(NA_real_)
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  min(sub)
}

#' For each i in sel, return its nearest-other-selected distance.
#' @keywords internal
.GprNearestInSel <- function(d, sel) {
  m <- length(sel)
  if (m < 2L) return(rep(Inf, m))
  sub <- d[sel, sel, drop = FALSE]
  diag(sub) <- Inf
  apply(sub, 1L, min)
}

#' Count pairs at the minimum distance (used by extended-improvement LS).
#' @keywords internal
.GprMinPairCount <- function(d, sel, dstar) {
  if (length(sel) < 2L) return(0L)
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
    cand_mask <- is.finite(g) & g > -Inf
    # candidates are points with g > -Inf (i.e., not yet selected)
    cand_idx <- which(g > -Inf)
    gv <- g[cand_idx]
    gmax <- max(gv)
    gmin <- min(gv)
    thresh <- gmin + alpha * (gmax - gmin)
    rcl <- cand_idx[gv >= thresh]
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
  if (m < 2L) return(sel)
  repeat {
    sel_sorted <- sort(sel)
    di <- .GprNearestInSel(d, sel_sorted)
    dstar <- min(di)
    pair_count <- .GprMinPairCount(d, sel_sorted, dstar)
    critical <- which(di <= dstar)  # indices into sel_sorted
    in_sel <- logical(n); in_sel[sel_sorted] <- TRUE
    out_idx <- which(!in_sel)
    best_swap <- NULL
    best_dstar <- dstar
    best_paircount <- pair_count
    for (ci in critical) {
      drop <- sel_sorted[ci]
      remaining <- sel_sorted[-ci]
      for (s in out_idx) {
        cand <- c(remaining, s)
        # quick: compute d* directly
        sub <- d[cand, cand, drop = FALSE]
        diag(sub) <- Inf
        new_dstar <- min(sub)
        if (new_dstar > best_dstar) {
          best_swap <- list(drop = drop, add = s)
          best_dstar <- new_dstar
          best_paircount <- sum(sub <= new_dstar) %/% 2L
        } else if (new_dstar == best_dstar) {
          new_paircount <- sum(sub <= new_dstar) %/% 2L
          if (new_paircount < best_paircount) {
            best_swap <- list(drop = drop, add = s)
            best_paircount <- new_paircount
          }
        }
      }
    }
    if (is.null(best_swap)) break
    sel <- c(sel_sorted[sel_sorted != best_swap$drop], best_swap$add)
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
  to_drop <- setdiff(x, y)  # |to_drop| = r
  to_add  <- setdiff(y, x)
  r <- length(to_drop)
  pk <- x
  best_sel <- x
  best_z <- .GprObjective(d, x)
  z_y <- .GprObjective(d, y)
  if (z_y > best_z) { best_sel <- y; best_z <- z_y }
  intermediates <- 0L
  for (k in seq_len(r)) {
    # pk_minus_y = to_drop intersect pk; pk_plus_y = to_add \ pk
    drop_cands <- intersect(pk, to_drop)
    add_cands  <- setdiff(to_add, pk)
    best_pair <- NULL
    best_pair_z <- -Inf
    # tie-breaker: lexicographic on (drop_index, add_index) -> iterate sorted
    for (i in sort(drop_cands)) {
      for (j in sort(add_cands)) {
        cand <- c(pk[pk != i], j)
        sub <- d[cand, cand, drop = FALSE]
        diag(sub) <- Inf
        zc <- min(sub)
        if (zc > best_pair_z) {
          best_pair_z <- zc
          best_pair <- c(i, j)
        }
      }
    }
    pk <- sort(c(pk[pk != best_pair[1L]], best_pair[2L]))
    intermediates <- intermediates + 1L
    if (best_pair_z > best_z) {
      best_sel <- pk
      best_z <- best_pair_z
    }
  }
  list(best = best_sel, objective = best_z, intermediates = intermediates)
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
.GprTryInsert <- function(d, ES, ES_z, sel, sel_z, dth) {
  z1 <- ES_z[1L]
  zb <- ES_z[length(ES_z)]
  hamm <- .GprHammingToES(sel, ES)
  dmin <- min(hamm)
  accept <- FALSE
  if (sel_z > z1) accept <- TRUE
  else if (sel_z > zb && dmin >= dth) accept <- TRUE
  # Don't admit a duplicate
  if (dmin == 0L) accept <- FALSE
  if (!accept) return(list(ES = ES, ES_z = ES_z, changed = FALSE))
  # remove the closest member (smallest Hamming distance); break ties on lowest z
  closest <- which(hamm == dmin)
  if (length(closest) > 1L) {
    closest <- closest[which.min(ES_z[closest])]
  } else {
    closest <- closest[1L]
  }
  ES <- ES[-closest]
  ES_z <- ES_z[-closest]
  # insert sorted (descending z)
  pos <- sum(ES_z >= sel_z) + 1L
  if (pos > length(ES)) {
    ES <- c(ES, list(sel))
    ES_z <- c(ES_z, sel_z)
  } else {
    ES <- c(ES[seq_len(pos - 1L)], list(sel), ES[pos:length(ES)])
    ES_z <- c(ES_z[seq_len(pos - 1L)], sel_z, ES_z[pos:length(ES_z)])
  }
  list(ES = ES, ES_z = ES_z, changed = TRUE)
}


# -- Outer loop: GRASP with PR static variant (Fig. 4) -------------------

#' GRASP with Path Relinking for the Max-Min Diversity Problem
#'
#' Solves the Max-Min Diversity Problem (discrete p-dispersion) with the
#' GRASP / path-relinking metaheuristic of Resende, Marti, Gallego & Duarte
#' (2010), static variant (their Fig. 4): a randomised-greedy construction
#' with extended-improvement local search builds and maintains an elite set,
#' followed by a single pass of path relinking over all elite pairs. On the
#' application benchmark this attains the highest \eqn{T_k} of the methods in
#' this package, at correspondingly higher cost.
#'
#' **Deterministic termination.** The refinement loop stops after
#' `max_no_improve` consecutive GRASP iterations that fail to improve the best
#' elite objective (rather than after a wall-clock budget). Given a fixed
#' `seed`, the entire run — construction RNG, iteration count, and result — is
#' therefore reproducible and machine-independent; `set.seed(seed)` controls
#' the same random stream the compiled kernel consumes via R's RNG. An
#' optional `time_budget_s` ceiling is available as a safety cap, but using a
#' finite value reintroduces machine-dependence and is off by default.
#'
#' This is a **dense-matrix-only** method: it materialises and repeatedly
#' subsets the full \eqn{n \times n} distance matrix, so it is suited to
#' instances small enough to hold that matrix. It offers no coordinate or
#' column-oracle path. For the matrix-free regime where the dense matrix is
#' infeasible, use [DropAddTSPoints()] or [Gonzalez()] (coordinate or
#' distance-column oracle path), whose
#' \eqn{T_k} lands within roughly a percent on the benchmark while scaling to
#' far larger instances.
#'
#' @param d Either a `dist` object or a square symmetric numeric matrix.
#' @param m Integer subset size, `2 <= m <= nrow(d)`.
#' @param max_no_improve Integer; stop after this many consecutive GRASP
#'   iterations without an improvement to the best elite objective. The
#'   primary, deterministic stopping criterion. Default 100.
#' @param max_iter Optional integer hard cap on GRASP refinement iterations
#'   (excluding the elite-set construction). `NULL` (default) leaves
#'   `max_no_improve` in sole control.
#' @param elite_size Size of the elite set |ES|. Default 10.
#' @param alpha RCL threshold; `alpha = 1` is pure greedy, `alpha = 0`
#'   uniform random. Default 0.8.
#' @param time_budget_s Optional wall-clock ceiling in seconds. Default `Inf`
#'   (no ceiling, fully reproducible). A finite value caps runtime but makes
#'   the result machine-dependent.
#' @param seed Optional integer; if supplied, `set.seed(seed)` is called at
#'   entry. `GraspPR` is genuinely stochastic (randomised construction and RCL
#'   sampling), so the seed governs the trajectory and the returned selection.
#' @return A list with elements
#'   \describe{
#'     \item{indices}{Integer vector of length `m`, 1-based, sorted ascending.}
#'     \item{objective}{Achieved MaxMin objective \eqn{T_k}.}
#'     \item{time_s}{Wall-clock seconds spent.}
#'     \item{iters}{Number of GRASP refinement iterations executed.}
#'     \item{pr_calls}{Number of path-relinking pair-applications run.}
#'   }
#' @references
#' Resende MGC, Marti R, Gallego M, Duarte A (2010). GRASP and path relinking
#' for the max-min diversity problem. \emph{Computers & Operations Research}
#' 37(3):498-508. \doi{10.1016/j.cor.2008.05.011}
#'
#' @seealso [DropAddTS()] and [DropAddTSPoints()] for scalable refinement;
#'   [ExactMaxMin()] for the proven optimum on small instances.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' res <- GraspPR(dist(pts), m = 5L, max_no_improve = 20L, elite_size = 4L,
#'                seed = 1L)
#' res$indices
#' @export
GraspPR <- function(d, m, max_no_improve = 100L, max_iter = NULL,
                    elite_size = 10L, alpha = 0.8, time_budget_s = Inf,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  m <- as.integer(m)
  stopifnot(m >= 2L, m <= n)
  elite_size <- as.integer(elite_size)
  stopifnot(elite_size >= 1L)
  max_no_improve <- as.integer(max_no_improve)
  stopifnot(max_no_improve >= 1L)
  if (is.null(max_iter)) {
    max_iter <- .Machine$integer.max
  } else {
    max_iter <- as.integer(max_iter)
    stopifnot(max_iter >= 0L)
  }
  if (!is.numeric(time_budget_s) || length(time_budget_s) != 1L ||
      is.na(time_budget_s) || time_budget_s <= 0) {
    stop("`time_budget_s` must be a single positive numeric (or Inf)")
  }

  out <- GraspPR_cpp(d, m, max_no_improve, max_iter, elite_size,
                     as.double(alpha), as.double(time_budget_s))
  list(
    indices   = sort(as.integer(out$indices)),
    objective = as.numeric(out$objective),
    time_s    = as.numeric(out$time_s),
    iters     = as.integer(out$iters),
    pr_calls  = as.integer(out$pr_calls)
  )
}

# Pure-R reference implementation of GraspPR, used as the parity oracle for
# the compiled kernel (see tests/testthat/test-gpr.R). It mirrors
# GraspPR_cpp() step for step, including the construction RNG draws, so that
# from a common `set.seed()` the two agree bit for bit. Not exported; callers
# use GraspPR().
#
# Termination is the deterministic stagnation rule: stop after
# `max_no_improve` consecutive iterations that do not raise the best elite
# objective ES_z[1] (which is monotone non-decreasing under .GprTryInsert()).
# `max_iter` is an optional hard cap; `time_budget_s` an optional ceiling
# (Inf = off) that leaves the result reproducible.
#' @keywords internal
.GraspPR_R <- function(d, m, max_no_improve, max_iter = .Machine$integer.max,
                       elite_size = 10L, alpha = 0.8, time_budget_s = Inf) {
  d <- .AsDistMatrix(d)
  n <- nrow(d)
  m <- as.integer(m)
  elite_size <- as.integer(elite_size)
  dth <- 5L
  t0 <- Sys.time()
  elapsed <- function() as.numeric(Sys.time() - t0, units = "secs")
  gated <- is.finite(time_budget_s)

  # Phase A: build initial elite set.
  ES <- vector("list", elite_size)
  ES_z <- numeric(elite_size)
  for (b in seq_len(elite_size)) {
    x  <- .GprConstruct(d, m, alpha)
    xp <- .GprLocalSearch(d, x)
    ES[[b]] <- xp
    ES_z[b] <- .GprObjective(d, xp)
  }
  ord <- order(ES_z, decreasing = TRUE)
  ES <- ES[ord]; ES_z <- ES_z[ord]

  # Phase B: GRASP iterations until `max_no_improve` consecutive non-improving
  # iterations (the deterministic criterion), an optional iteration cap, or an
  # optional wall-clock ceiling.
  iters <- 0L
  no_improve <- 0L
  best_z_B <- ES_z[1L]
  repeat {
    if (no_improve >= max_no_improve) break
    if (iters >= max_iter) break
    if (gated && elapsed() >= time_budget_s) break
    x  <- .GprConstruct(d, m, alpha)
    xp <- .GprLocalSearch(d, x)
    zp <- .GprObjective(d, xp)
    res <- .GprTryInsert(d, ES, ES_z, xp, zp, dth)
    ES <- res$ES; ES_z <- res$ES_z
    iters <- iters + 1L
    if (ES_z[1L] > best_z_B) {
      best_z_B <- ES_z[1L]
      no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
    }
  }

  # Phase C: path relinking over all elite pairs (deterministic; no RNG).
  best_sel <- ES[[1L]]
  best_z   <- ES_z[1L]
  pr_calls <- 0L
  k <- length(ES)
  if (k >= 2L && !(gated && elapsed() >= time_budget_s)) {
    done <- FALSE
    for (i in seq_len(k - 1L)) {
      if (done) break
      for (j in (i + 1L):k) {
        pr1 <- .GprPathRelink(d, ES[[i]], ES[[j]])
        pr2 <- .GprPathRelink(d, ES[[j]], ES[[i]])
        pr_calls <- pr_calls + 2L
        y_sel <- if (pr1$objective >= pr2$objective) pr1$best else pr2$best
        yp <- .GprLocalSearch(d, y_sel)
        zp <- .GprObjective(d, yp)
        if (zp > best_z) {
          best_z <- zp
          best_sel <- yp
        }
        if (gated && elapsed() >= time_budget_s) {
          done <- TRUE
          break
        }
      }
    }
  }

  list(
    indices   = sort(as.integer(best_sel)),
    objective = best_z,
    time_s    = elapsed(),
    iters     = iters,
    pr_calls  = pr_calls
  )
}
