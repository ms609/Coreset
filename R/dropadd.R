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
  row_sums <- rowSums(dmat)
  S[1L] <- which.max(row_sums)  # ties: which.max → smallest index (deterministic)
  in_S <- logical(n)
  in_S[S[1L]] <- TRUE

  # Streamlined records.
  # min_dist(x):       min over y in current_S of d(x, y), excluding x itself
  # sum_dist(x):       sum over y in current_S of d(x, y) (d(x,x)=0 so harmless)
  # min_dist_count(x): | { y in current_S, y != x : d(x, y) == min_dist(x) } |
  min_dist       <- dmat[, S[1L]]            # length n; min_dist(S[1])=0 (self)
  min_dist[S[1L]] <- Inf                     # never want to compare a selected
                                             # point to itself in min calcs
  sum_dist       <- dmat[, S[1L]]
  min_dist_count <- rep(1L, n)
  min_dist_count[S[1L]] <- 0L                # no "other selected" yet

  if (m >= 2L) {
    for (h in 2L:m) {
      # Among AddX = Z \ S, pick x maximising (min_dist, sum_dist) lex,
      # ties → smallest index.
      cand <- which(!in_S)
      md <- min_dist[cand]
      best_min <- max(md)
      tied <- cand[md == best_min]
      if (length(tied) > 1L) {
        sd <- sum_dist[tied]
        best_sum <- max(sd)
        tied <- tied[sd == best_sum]
      }
      x_new <- tied[1L]   # smallest index among lex-ties
      S[h] <- x_new
      in_S[x_new] <- TRUE

      # Update streamlined records for ADD (Algorithm 3).
      # For every x != x_new:
      d_col <- dmat[, x_new]
      others <- seq_len(n)[-x_new]
      d_to_new <- d_col[others]
      # sum_dist
      sum_dist[others] <- sum_dist[others] + d_to_new
      # Cases
      md_o <- min_dist[others]
      caseB <- d_to_new < md_o
      caseA <- d_to_new == md_o & !caseB
      # Apply
      if (any(caseB)) {
        idx <- others[caseB]
        min_dist[idx] <- d_to_new[caseB]
        min_dist_count[idx] <- 1L
      }
      if (any(caseA)) {
        idx <- others[caseA]
        min_dist_count[idx] <- min_dist_count[idx] + 1L
      }
      # x_new itself: its min_dist becomes min over (S - {x_new}) of d.
      # Since we'd already maintained min_dist[x_new] = Inf, we need to set
      # it now to the min distance from x_new to the previously-selected
      # points (which is what its min_dist would have been without the
      # Inf override).
      prev_S <- S[seq_len(h - 1L)]
      d_xnew <- dmat[x_new, prev_S]
      mn <- min(d_xnew)
      min_dist[x_new] <- mn
      min_dist_count[x_new] <- sum(d_xnew == mn)
    }
  }

  list(S = S, in_S = in_S, min_dist = min_dist, sum_dist = sum_dist,
       min_dist_count = min_dist_count)
}

# nocov start
# MMDPo objective from streamlined records: for the current set X,
#   maxmin = min over x in X of min_dist(x)
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

# Recompute min_dist(x) and min_dist_count(x) from scratch given current S
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
#' \eqn{(\mathrm{min\_dist}, \mathrm{sum\_dist})} is added. The FIFO
#' invariant guarantees that across any window of \eqn{m} iterations every
#' initially-selected point is dropped exactly once before any re-eviction.
#'
#' @param d A \code{dist} object or square symmetric numeric matrix.
#' @param m Integer; subset size, \eqn{2 \le m \le n}.
#' @param max_no_improve Integer; stop after this many consecutive drop-add
#'   iterations that do not improve the best objective. The primary,
#'   deterministic stopping criterion. The search is RNG-free (ties broken by
#'   smallest index), so for a given instance the result is reproducible and
#'   machine-independent. Default 5000.
#' @param max_iter Optional integer hard cap on iterations (excluding
#'   construction). \code{NULL} (default) leaves \code{max_no_improve} in sole
#'   control.
#' @param time_budget_s Optional wall-clock ceiling in seconds, checked at
#'   iteration boundaries. Default \code{Inf} (no ceiling, fully reproducible).
#'   A finite value caps runtime but makes the result machine-dependent.
#' @param seed Optional integer; if non-NULL, \code{set.seed(seed)} is called
#'   at entry. The algorithm is deterministic up to ties (broken by smallest
#'   index), so the seed has no observable effect on the solution; it is
#'   exposed for API parity with stochastic methods.
#' @param progress Logical; show a start/done status line. Default: `TRUE` in
#'   interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`). No effect when
#'   `.verify = TRUE` (testing path).
#' @param .verify Logical (testing only); if `TRUE`, routes to the R reference
#'   loop and brute-force asserts the streamlined records at every iteration.
#'   Default `FALSE` (the C++ fast path).
#' @param .trace Optional environment (testing only); if supplied, the dropped
#'   and added index sequences are written into it as `drops` and `adds`.
#'
#' @return List with elements
#'   \describe{
#'     \item{indices}{integer(m), 1-based selected indices sorted ascending.}
#'     \item{objective}{numeric(1), achieved MaxMin objective
#'       \eqn{\min_{i \ne j \in S} d_{ij}}.}
#'     \item{secondary}{numeric(1), achieved sum of pairwise distances over
#'       \eqn{S} (upper-triangle sum).}
#'     \item{time_s}{numeric(1), wall-clock seconds spent.}
#'     \item{iters}{integer(1), main-loop iterations executed (excluding the
#'       construction phase).}
#'   }
#'
#' @references \insertAllCited{}
#'
#' @export
DropAddTS <- function(d, m, max_no_improve = 5000L, max_iter = NULL,
                      time_budget_s = Inf, seed = NULL,
                      progress = getOption("MaxMin.progress", interactive()),
                      .verify = FALSE, .trace = NULL) {
  # .verify: if TRUE, brute-force recompute and assert all streamlined
  #   records (min_dist, sum_dist, min_dist_count) at every iteration.
  #   Routes to the R reference loop; off by default.
  # .trace: optional environment; if supplied, the algorithm writes
  #   `drops` (integer vector of dropped indices in order) and `adds`
  #   (integer vector of added indices in order) into it. Used by the
  #   FIFO test to inspect drop order without re-implementing the loop.
  #   Supported by both the R and C++ paths.
  if (!is.null(seed)) set.seed(seed)
  dmat <- .AsDistMatrix(d)
  n <- nrow(dmat)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L || m > n) {
    stop("`m` must be a single integer with 2 <= m <= nrow(d)")
  }
  max_no_improve <- as.integer(max_no_improve)
  if (length(max_no_improve) != 1L || is.na(max_no_improve) ||
      max_no_improve < 1L) {
    stop("`max_no_improve` must be a single positive integer")
  }
  if (!is.null(max_iter)) {
    max_iter <- as.integer(max_iter)
    if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 0L) {
      stop("`max_iter` must be NULL or a single non-negative integer")
    }
  }
  if (!is.numeric(time_budget_s) || length(time_budget_s) != 1L ||
      is.na(time_budget_s) || time_budget_s <= 0) {
    stop("`time_budget_s` must be a single positive numeric (or Inf)")
  }

  t0 <- Sys.time()
  eps <- 1e-9

  # --- C++ fast path. .verify routes to R for the brute-force assertion. --
  if (!.verify) {
    cpp_max_iter <- if (is.null(max_iter)) .Machine$integer.max else max_iter
    want_trace <- !is.null(.trace)
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (n = {n}, m = {m}, budget = {time_budget_s}s)",
        .auto_close = FALSE
      )
    }
    out <- DropAddTS_cpp(dmat, m, as.double(time_budget_s),
                         cpp_max_iter, max_no_improve, want_trace)
    if (want_trace) {
      .trace$drops <- out$drops
      .trace$adds  <- out$adds
    }
    time_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (progress) {
      iters_msg <- as.integer(out$iters)
      tk_msg    <- as.numeric(out$objective)
      cli::cli_process_done(
        msg = "DropAdd: {iters_msg} iters, T_k = {signif(tk_msg, 4)}, {round(time_s, 1)}s"
      )
    }
    return(list(
      indices   = sort(as.integer(out$indices)),
      objective = as.numeric(out$objective),
      secondary = as.numeric(out$secondary),
      time_s    = time_s,
      iters     = as.integer(out$iters)
    ))
  }

  # --- Construction (Algorithm 1) -----------------------------------------
  cons <- .DropAddConstruct(dmat, m)
  S        <- cons$S
  in_S     <- cons$in_S
  min_dist <- cons$min_dist
  sum_dist <- cons$sum_dist
  min_dist_count <- cons$min_dist_count

  # Track best-known MMDPo solution. Compute the initial objective from the
  # streamlined records (not via .MMDPoObj's upper-triangle sum) so that the
  # R fallback path and the C++ port are bit-identical from the first iter:
  # `.MMDPoObj` sums `dmat[S, S]` upper-tri in row-major order, the
  # streamlined `sum(sum_dist[S]) / 2` accumulates m partial sums — both
  # compute the same quantity but differ in the last ULP.
  best_S <- S
  best_maxmin  <- min(min_dist[S])
  best_sumpair <- sum(sum_dist[S]) / 2
  best_score   <- best_maxmin + eps * best_sumpair

  # --- Drop-Add tabu search (Algorithm 2) ---------------------------------
  # FIFO via circular buffer: S is a length-m vector treated as a queue.
  # `head` indexes the oldest member. Each iteration drops S[head], writes
  # the new x* into S[head], and advances head modulo m. The Porumbel
  # iter_stamp / which.min(iter_stamp[S]) pair is redundant under this
  # invariant — head IS the FIFO head.
  head        <- 1L
  iters_done  <- 0L
  no_improve  <- 0L
  t0_num      <- unclass(t0)
  check_every <- 64L
  countdown   <- check_every
  effective_max <- if (is.null(max_iter)) .Machine$integer.max else max_iter
  if (m >= n) effective_max <- 0L  # all points selected: no drop-add move exists
  if (!is.null(.trace)) {
    .trace$drops <- integer(0)
    .trace$adds  <- integer(0)
  }

  # Internal brute-force assertion of streamlined records.
  .verify_records <- function(S, min_dist, min_dist_count, sum_dist) {
    for (xx in seq_len(n)) {
      others <- if (xx %in% S) setdiff(S, xx) else S
      if (!length(others)) next  # nocov
      vals <- dmat[xx, others]
      tm <- min(vals); tc <- sum(vals == tm); ts <- sum(dmat[xx, S])
      stopifnot(abs(min_dist[xx] - tm) < 1e-9,
                min_dist_count[xx] == tc,
                abs(sum_dist[xx] - ts) < 1e-6)
    }
    invisible(TRUE)
  }
  if (.verify) .verify_records(S, min_dist, min_dist_count, sum_dist)

  repeat {
    # Termination checks. max_iter is cheap (integer compare); Sys.time is
    # throttled to once per `check_every` iterations because POSIXct
    # creation + difftime dispatch was ~15% of per-iter cost.
    if (iters_done >= effective_max) break
    if (no_improve >= max_no_improve) break
    countdown <- countdown - 1L
    if (countdown == 0L) {
      if (unclass(Sys.time()) - t0_num >= time_budget_s) break  # nocov
      countdown <- check_every
    }

    # 1. DROP: head IS the FIFO head under the circular-buffer invariant.
    x_hash <- S[head]

    # Apply drop to streamlined records (Algorithm 4). Full-vector updates
    # operate on all n entries; dmat[x_hash, x_hash] = 0 so sum_dist[x_hash]
    # is unchanged, and x_hash's min_dist/min_dist_count get clobbered by
    # the trailing self-recompute regardless.
    in_S[x_hash] <- FALSE
    d_drop_full <- dmat[, x_hash]
    sum_dist <- sum_dist - d_drop_full
    affected_mask <- d_drop_full <= min_dist
    affected_mask[x_hash] <- FALSE         # x_hash gets its own recompute below
    if (any(affected_mask)) {
      min_dist_count[affected_mask] <- min_dist_count[affected_mask] - 1L
    }
    need_recompute <- which(affected_mask & min_dist_count == 0L)
    S_after_drop <- S[-head]
    if (length(need_recompute)) {
      # Vectorised recompute. `sub` is k × (m-1); rows for need_recompute
      # entries still in S need their self entry masked to Inf.
      sub <- dmat[need_recompute, S_after_drop, drop = FALSE]
      in_S_rows <- in_S[need_recompute]
      if (any(in_S_rows)) {
        sel <- which(in_S_rows)
        cols <- match(need_recompute[sel], S_after_drop)
        sub[cbind(sel, cols)] <- Inf
      }
      mns  <- as.numeric(do.call(pmin.int, asplit(sub, 2L)))
      cnts <- rowSums(sub == mns)
      finite_mns <- is.finite(mns)
      if (all(finite_mns)) {
        min_dist[need_recompute]       <- mns
        min_dist_count[need_recompute] <- cnts
      } else {
        bad <- !finite_mns
        min_dist[need_recompute[bad]]       <- Inf
        min_dist_count[need_recompute[bad]] <- 0L
        if (any(finite_mns)) {
          good <- need_recompute[finite_mns]
          min_dist[good]       <- mns[finite_mns]
          min_dist_count[good] <- cnts[finite_mns]
        }
      }
    }
    # x_hash's own min_dist record: distance to nearest surviving peer.
    if (length(S_after_drop)) {
      row_vals <- dmat[x_hash, S_after_drop]
      mn <- min(row_vals)
      min_dist[x_hash] <- mn
      min_dist_count[x_hash] <- sum(row_vals == mn)
    } else {
      min_dist[x_hash] <- Inf         # nocov
      min_dist_count[x_hash] <- 0L   # nocov
    }

    # 2. ADD: argmax over Add X(k) = Z - X(k) of (min_dist, sum_dist), ties →
    # smallest idx. x_hash is excluded for this iteration (Porumbel et al. 2011,
    # p.281): the just-dropped point cannot be re-added immediately — the tabu
    # rule that prevents looping. It is eligible again once head advances.
    cand <- which(!in_S)
    cand <- cand[cand != x_hash]
    md <- min_dist[cand]
    best_min <- max(md)
    tied <- cand[md == best_min]
    if (length(tied) > 1L) {
      sd <- sum_dist[tied]
      best_sum <- max(sd)
      tied <- tied[sd == best_sum]
    }
    x_new <- tied[1L]

    # Update streamlined records for ADD (Algorithm 3). Same full-vector
    # discipline; x_new's own record is overwritten by the trailing block.
    in_S[x_new] <- TRUE
    d_col_full <- dmat[, x_new]
    sum_dist <- sum_dist + d_col_full
    caseB <- d_col_full <  min_dist
    caseA <- d_col_full == min_dist        # < and == are disjoint; no !caseB mask needed
    if (any(caseB)) {
      min_dist[caseB]       <- d_col_full[caseB]
      min_dist_count[caseB] <- 1L
    }
    if (any(caseA)) {
      min_dist_count[caseA] <- min_dist_count[caseA] + 1L
    }
    # x_new's own min_dist: distance to peers in S \ {x_hash}.
    row_vals2 <- dmat[x_new, S_after_drop]
    mn2 <- min(row_vals2)
    min_dist[x_new] <- mn2
    min_dist_count[x_new] <- sum(row_vals2 == mn2)

    # Write x_new into the head slot and advance.
    S[head] <- x_new
    head <- if (head == m) 1L else head + 1L

    # 3. Test for improvement of best-known MMDPo solution.
    cur_maxmin  <- min(min_dist[S])
    cur_sumpair <- sum(sum_dist[S]) / 2     # double-counted
    cur_score   <- cur_maxmin + eps * cur_sumpair
    if (cur_score > best_score) {
      best_S       <- S
      best_maxmin  <- cur_maxmin
      best_sumpair <- cur_sumpair
      best_score   <- cur_score
      no_improve   <- 0L
    } else {
      no_improve   <- no_improve + 1L
    }

    iters_done <- iters_done + 1L

    if (!is.null(.trace)) {
      .trace$drops <- c(.trace$drops, x_hash)
      .trace$adds  <- c(.trace$adds, x_new)
    }
    if (.verify) .verify_records(S, min_dist, min_dist_count, sum_dist)
  }

  time_s <- unclass(Sys.time()) - t0_num
  list(
    indices   = sort(as.integer(best_S)),
    objective = as.numeric(best_maxmin),
    secondary = as.numeric(best_sumpair),
    time_s    = time_s,
    iters     = as.integer(iters_done)
  )
}
