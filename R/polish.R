# polish.R

#' Local-search polish for a max-min diversity selection
#'
#' Given a selection `idx` of size `k >= 2`, search for 1-swaps that increase
#' \eqn{T_k = \min_{a \neq b \in S} d(a, b)} (or, at fixed `T_k`, reduce the
#' number of pairs achieving it). The neighbourhood is critical-edge-anchored:
#' the candidate removal `x` ranges over endpoints of pairs achieving the
#' current minimum, and the candidate insertion `w` over the `limit` nearest
#' neighbours of `x` in the full data (excluding the current selection).
#'
#' Acceptance follows Della Croce et al.'s flat-landscape rule: accept a swap
#' if `new_T > T` *or* `new_T == T` *and* the count of pairs at the minimum
#' strictly decreases. Both `T_k` and `-n_critical` are monotone non-decreasing
#' across accepted swaps, so iteration terminates without a tabu list.
#'
#' The intended use is as a post-step on any greedy max-min selector — for
#' example [Gonzalez()]. Calling `PolishSelection()` on a 1-swap local optimum
#' is a no-op.
#'
#' @param d A `dist` object or square symmetric numeric matrix of pairwise
#'   distances.
#' @param idx Integer vector of selected row/col indices (`k = length(idx)`).
#' @param limit Integer: maximum neighbour rank to scan per critical endpoint.
#'   Default `20L`. Smaller values are faster but may miss improving swaps.
#' @param max_passes Integer: hard cap on outer iterations (safety guard;
#'   theoretical termination is guaranteed by the lex-monotone objective).
#'   Default `200L`.
#' @param progress Logical; show a start/done status line. Default: `TRUE` in
#'   interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
#' @return Integer vector of the same length as `idx`, with attributes
#'   `"passes"` and `"swaps"` (both equal to the number of accepted swaps).
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d   <- dist(pts)
#' s0  <- Gonzalez(d, 5L, first = 1L)
#' s1  <- PolishSelection(d, s0)
#' TkScore(d, s1) >= TkScore(d, s0)
#' @export
PolishSelection <- function(d, idx, limit = 20L, max_passes = 200L,
                            progress = getOption("MaxMin.progress", interactive())) {
  d   <- .AsDistMatrix(d)
  idx <- as.integer(idx)
  if (length(idx) < 2L) return(idx)
  if (progress) {
    cli::cli_process_start(
      "Polishing selection (k = {length(idx)})",
      .auto_close = FALSE
    )
  }
  t0  <- Sys.time()
  out <- PolishMaximin_cpp(d, idx, as.integer(limit), as.integer(max_passes))
  if (progress) {
    time_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .passes <- attr(out, "passes")
    .swaps  <- attr(out, "swaps")
    cli::cli_process_done(
      msg = "Polish: {.passes} passes, {.swaps} swaps, {round(time_s, 2)}s"
    )
  }
  out
}
