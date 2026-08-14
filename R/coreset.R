# coreset.R
#
# Composable-coreset candidate thinning for the heavy MMDP solvers.
#
# The heavy solvers (DropAdd tabu search, GRASP with path-relinking) do not
# scale to large N: DropAdd's wall-time is prohibitive and the matrix-only
# GRASP needs a dense N x N matrix. The standard fix is a composable coreset
# (Indyk et al. 2014; Aghamolaei et al. 2015): run the cheap farthest-first
# greedy (FarFirst) to reduce the N candidates to a small intermediate subset of
# m points (m << N, m >= k), then run the expensive solver on just those m and
# map the chosen indices back to the original numbering.
#
# This file holds the two solver-agnostic helpers behind the `maxCandidates`
# parameter of DropAdd() and Grasp(): `.ResolveCap()` (validate/normalise the
# cap) and `.FarFirstThin()` (build the coreset, run the solver on it, and
# perform the index round-trip). The per-solver wiring lives in each solver.

#' Validate and normalise a `maxCandidates` thinning cap
#'
#' Shared by [DropAdd()] and [Grasp()]. Decides whether candidate thinning is
#' active and, if so, the intermediate coreset size `m` to thin to.
#'
#' @param maxCandidates The user-supplied cap: a positive integer (thin to this
#'   many candidates when it is below `n`), or `0` / `Inf` to disable thinning.
#' @param n Integer: the number of candidate points in the full problem.
#' @param k Integer: the target subset size.
#' @return `.ResolveCap()` returns `NA_integer_` when thinning is disabled
#'   (`0` / `Inf`) or non-binding (`maxCandidates >= n`); otherwise the integer
#'   coreset size `m` (`k <= m < n`). Errors on a non-integer, negative, `NA`,
#'   or non-scalar cap, or a positive cap below `k`.
#' @keywords internal
.ResolveCap <- function(maxCandidates, n, k) {
  if (length(maxCandidates) != 1L || is.na(maxCandidates)) {
    stop("`maxCandidates` must be a single number (0 or Inf disables thinning)")
  }
  # 0 / Inf disable thinning. (Test the double before as.integer() so a "huge
  # off-switch" value cannot overflow the integer range.)
  if (is.infinite(maxCandidates) || maxCandidates == 0) {
    return(NA_integer_)
  }
  if (maxCandidates < 0 || maxCandidates != floor(maxCandidates)) {
    stop("`maxCandidates` must be a positive integer, 0, or Inf")
  }
  # Cap at or above n is not binding: run on the full problem (a no-op coreset).
  if (maxCandidates >= n) {
    return(NA_integer_)
  }
  m <- as.integer(maxCandidates)
  if (m < k) {
    stop("`maxCandidates` (", m, ") must be >= k (", k, ")")
  }
  m
}

#' Run a solver on a farthest-first coreset and map indices back
#'
#' Implements the composable-coreset path of [DropAdd()] / [Grasp()] (dispatched
#' there when their `maxCandidates` cap binds). It builds an `m`-point coreset
#' with [FarFirst()], restricts the problem to those `m` points, runs the
#' supplied solver on the restriction, and maps the returned indices back to the
#' original numbering.
#'
#' Only the integer index values returned by the solver need remapping.
#' The result is sorted ascending.
#'
#' @param k Integer: target subset size.
#' @param m Integer: coreset size (`k <= m < N`, as returned by [.ResolveCap()]).
#' @param d Square distance matrix of the full problem, or `NULL` on the
#'   coordinate path.
#' @param points `N x dim` coordinate matrix of the full problem, or `NULL` on
#'   the distance-matrix path. Exactly one of `d` / `points` is non-`NULL`.
#' @param RunOnSubset Function `function(d, points)` that runs the downstream
#'   solver on the restricted problem and returns its `MaxMinSelection`. It is
#'   called with the restricted distance matrix (`d = d[core, core]`) on the
#'   matrix path, or the restricted coordinates (`points = points[core, ]`) on
#'   the coordinate path; the unused argument is `NULL`.
#' @param label Character naming the calling solver, used in the thinning
#'   warning (e.g. `"DropAdd"`).
#' @return `.FarFirstThin()` returns the solver's `MaxMinSelection` with its
#'   indices mapped to original-space row indices (sorted ascending).
#' @keywords internal
.FarFirstThin <- function(k, m, d = NULL, points = NULL, RunOnSubset, label) {
  N <- if (is.null(points)) nrow(d) else nrow(points)
  # Deterministic peripheral seed: a single bare Gonzalez pass, no RNG drawn.
  core <- as.integer(FarFirst(m, d = d, points = points, strategy = "peripheral"))
  warning(label, ": farthest-first thinned ", N, " candidates to ", m,
          " (maxCandidates); pass maxCandidates = 0 to disable", call. = FALSE)
  local <- if (is.null(points)) {
    RunOnSubset(d = d[core, core, drop = FALSE], points = NULL)
  } else {
    RunOnSubset(d = NULL, points = points[core, , drop = FALSE])
  }
  # Round-trip: map coreset-local index values to original-space indices.
  mapped <- core[as.integer(local)]
  out <- mapped[order(mapped)]
  attributes(out) <- attributes(local)
  out
}
