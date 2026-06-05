# dropadd_points.R
#
# Matrix-free (coordinate-based) DropAdd Tabu Search for the Max-Min Diversity
# Problem (Porumbel, Hao & Glover 2011).
#
# Coordinate counterpart of DropAddTS() in competitors_dropadd.R. The matrix
# wrapper coerces its input to a dense n x n distance matrix and feeds it to
# DropAddTS_cpp; this wrapper instead passes the raw n x dim coordinate matrix
# to DropAddTS_points_cpp, which recomputes each needed distance column on the
# fly. The dense matrix is never built, so the SOTA heuristic runs at n far
# beyond the matrix path's ceiling (R's as.matrix.dist overflows at n = 46340,
# and an n = 58000 double matrix is ~27 GB).
#
# The algorithm is identical to DropAddTS() in every respect except the seed:
# Porumbel's argmax-row-sum seed is O(n^2 * dim), so the kernel substitutes the
# O(n * dim) "farthest point from the centroid" proxy (documented in
# src/dropadd_mf.cpp). On instances where the two seed rules coincide, the
# entire trajectory matches the matrix path; otherwise quality is comparable.

#' Matrix-free DropAdd Tabu Search for the Max-Min Diversity Problem
#'
#' Coordinate-based counterpart of [DropAddTS()] that never materialises the
#' dense \eqn{n \times n} distance matrix. Each needed distance column
#' \eqn{d(\cdot, x)} is recomputed from the supplied coordinates on the fly in
#' \eqn{O(n \cdot \mathrm{dim})}, giving \eqn{O(n)} working memory. This lets the
#' SOTA MaxMin heuristic of Porumbel, Hao & Glover (2011) run on point sets far
#' larger than the matrix path can hold (R's `as.matrix.dist` overflows at
#' \eqn{n = 46340}; an \eqn{n = 58000} double matrix is roughly 27 GB).
#'
#' The construction (Algorithm 1), FIFO drop-add tabu search (Algorithm 2), and
#' streamlined neighbour evaluation (Algorithms 3-4) are identical to
#' [DropAddTS()], including the exclusion of the just-dropped point from the add
#' candidates for one iteration (Porumbel et al. 2011, p.281). The MMDPo
#' objective optimised is
#' \deqn{\min_{x,y \in X} d(x,y) + \epsilon \sum_{x,y \in X} d(x,y),}
#' with \eqn{\epsilon = 10^{-9}}.
#'
#' Seed deviation. Porumbel seeds the construction with the maximum-row-sum
#' point, \eqn{\arg\max_x \sum_y d(x,y)}, which costs \eqn{O(n^2 \cdot
#' \mathrm{dim})} — the very expense this variant avoids. We substitute the
#' \eqn{O(n \cdot \mathrm{dim})} proxy \eqn{\arg\max_x \lVert x - \bar{x}
#' \rVert} (the point farthest from the coordinate centroid), which
#' approximates the peripheral max-row-sum point. The remainder of the
#' algorithm is faithful, so on instances where the two seed rules coincide the
#' trajectory matches [DropAddTS()] exactly; otherwise the final quality is
#' comparable (within a few percent on the MaxMin objective, empirically).
#'
#' Distances reproduce [stats::dist()]'s Euclidean bits exactly (see
#' `src/maximin_points.cpp` for the bit-equivalence argument), so on data where
#' the seeds coincide the matrix-free and matrix paths are bit-identical under a
#' toolchain whose `stats` package and this kernel use the same floating-point
#' contraction settings.
#'
#' @param points A numeric \eqn{n \times \mathrm{dim}} coordinate matrix (or an
#'   object coercible to one via `as.matrix`). Must be complete (no `NA`); the
#'   coordinate path is defined only for complete Euclidean data.
#' @param m Integer; subset size, \eqn{2 \le m \le n}.
#' @param max_no_improve Integer; stop after this many consecutive drop-add
#'   iterations that do not improve the best objective. The primary,
#'   deterministic stopping criterion. The search is RNG-free (ties broken by
#'   smallest index), so the result is reproducible and machine-independent.
#'   Default 5000.
#' @param max_iter Optional integer hard cap on main-loop iterations (excluding
#'   construction). `NULL` (default) leaves `max_no_improve` in sole control.
#' @param time_budget_s Optional wall-clock ceiling in seconds, checked
#'   periodically at iteration boundaries. Default `Inf` (no ceiling, fully
#'   reproducible). A finite value caps runtime but makes the result
#'   machine-dependent.
#' @param seed Optional integer; if non-`NULL`, `set.seed(seed)` is called at
#'   entry. The algorithm is deterministic up to ties (broken by smallest
#'   index), so the seed has no observable effect on the solution; it is exposed
#'   for API parity with stochastic methods and with [DropAddTS()].
#' @param progress Logical; show a start/done status line. Default: `TRUE` in
#'   interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
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
#' @references
#' Porumbel D, Hao J-K, Glover F (2011). A simple and effective algorithm for
#' the MaxMin diversity problem. \emph{Annals of Operations Research}
#' 186:275-293.
#'
#' @seealso [DropAddTS()] for the matrix-based path used on smaller instances.
#' @export
DropAddTSPoints <- function(points, m, max_no_improve = 5000L, max_iter = NULL,
                            time_budget_s = Inf, seed = NULL,
                            progress = getOption("MaxMin.progress", interactive())) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  points <- .AsPointsMatrix(points)
  n <- nrow(points)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L || m > n) {
    stop("`m` must be a single integer with 2 <= m <= nrow(points)")
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
  cpp_max_iter <- if (is.null(max_iter)) .Machine$integer.max else max_iter
  if (progress) {
    cli::cli_process_start(
      "DropAdd tabu search (n = {n}, m = {m}, budget = {time_budget_s}s)",
      .auto_close = FALSE
    )
  }
  out <- DropAddTS_points_cpp(points, m, as.double(time_budget_s),
                              cpp_max_iter, max_no_improve, FALSE)
  time_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (progress) {
    .iters <- as.integer(out$iters)
    .tk    <- as.numeric(out$objective)
    cli::cli_process_done(
      msg = "DropAdd: {.iters} iters, T_k = {signif(.tk, 4)}, {round(time_s, 1)}s"
    )
  }

  # Return:
  list(
    indices   = sort(as.integer(out$indices)),
    objective = as.numeric(out$objective),
    secondary = as.numeric(out$secondary),
    time_s    = time_s,
    iters     = as.integer(out$iters)
  )
}
