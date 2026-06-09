# dropadd_points.R
#
# Matrix-free (coordinate-based) DropAdd Tabu Search for the Max-Min Diversity
# Problem (Porumbel, Hao & Glover 2011).
#
# Coordinate counterpart of DropAdd() in competitors_dropadd.R. The matrix
# wrapper coerces its input to a dense n x n distance matrix and feeds it to
# DropAdd_cpp; this wrapper instead passes the raw n x dim coordinate matrix
# to DropAdd_points_cpp, which recomputes each needed distance column on the
# fly. The dense matrix is never built, so the SOTA heuristic runs at n far
# beyond the matrix path's ceiling (R's as.matrix.dist overflows at n = 46340,
# and an n = 58000 double matrix is ~27 GB).
#
# The algorithm is identical to DropAdd() in every respect except the seed:
# Porumbel's argmax-row-sum seed is O(n^2 * dim), so the kernel substitutes the
# O(n * dim) "farthest point from the centroid" proxy (documented in
# src/dropadd_mf.cpp). On instances where the two seed rules coincide, the
# entire trajectory matches the matrix path; otherwise quality is comparable.

#' Matrix-free DropAdd Tabu Search for the Max-Min Diversity Problem
#'
#' Coordinate-based counterpart of [DropAdd()] that never materialises the
#' dense \eqn{n \times n} distance matrix. Each needed distance column
#' \eqn{d(\cdot, x)} is recomputed from the supplied coordinates on the fly in
#' \eqn{O(n \cdot \mathrm{dim})}, giving \eqn{O(n)} working memory. This lets the
#' SOTA MaxMin heuristic of \insertCite{Porumbel2011;textual}{MaxMin} run on point sets far
#' larger than the matrix path can hold (R's `as.matrix.dist` overflows at
#' \eqn{n = 46340}; an \eqn{n = 58000} double matrix is roughly 27 GB).
#'
#' The construction (Algorithm 1), FIFO drop-add tabu search (Algorithm 2), and
#' streamlined neighbour evaluation (Algorithms 3-4) are identical to
#' [DropAdd()], including the exclusion of the just-dropped point from the add
#' candidates for one iteration (\insertCite{Porumbel2011;textual}{MaxMin}, p.281). The MMDPo
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
#' trajectory matches [DropAdd()] exactly; otherwise the final quality is
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
#' @param plateau Integer; stop after this many consecutive drop-add
#'   iterations that do not improve the best objective. The primary,
#'   deterministic stopping criterion. The search is RNG-free (ties broken by
#'   smallest index), so the result is reproducible and machine-independent.
#'   Default 5000.
#' @param maxIter Optional integer hard cap on main-loop iterations (excluding
#'   construction). `NULL` (default) leaves `plateau` in sole control.
#' @param timeBudgetS Optional wall-clock ceiling in seconds, checked
#'   periodically at iteration boundaries. Default `Inf` (no ceiling, fully
#'   reproducible). A finite value caps runtime but makes the result
#'   machine-dependent.
#' @param seed Optional integer; if non-`NULL`, `set.seed(seed)` is called at
#'   entry. The algorithm is deterministic up to ties (broken by smallest
#'   index), so the seed has no observable effect on the solution; it is exposed
#'   for API parity with stochastic methods and with [DropAdd()].
#' @param progress Logical; show a start/done status line. Default: `TRUE` in
#'   interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
#'
#' @return An integer vector of length \code{m} containing the 1-based selected
#'   indices (sorted ascending), with attributes:
#'   \describe{
#'     \item{score}{numeric(1), achieved MaxMin objective
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
#' @seealso [DropAdd()] for the matrix-based path used on smaller instances.
#' @export
DropAddPoints <- function(points, m, plateau = 5000L, maxIter = NULL,
                          timeBudgetS = Inf, seed = NULL,
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

  t0 <- Sys.time()
  cppMaxIter <- if (is.null(maxIter)) .Machine$integer.max else maxIter
  if (progress) {
    cli::cli_process_start(
      "DropAdd tabu search (n = {n}, m = {m}, budget = {timeBudgetS}s)",
      .auto_close = FALSE
    )
  }
  out <- DropAdd_points_cpp(points, m, as.double(timeBudgetS),
                              cppMaxIter, plateau, FALSE)
  timeS <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (progress) {
    itersMsg <- as.integer(out$iters)
    tkMsg    <- as.numeric(out$score)
    cli::cli_process_done(
      msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
    )
  }

  structure(
    sort(as.integer(out$indices)),
    score     = as.numeric(out$score),
    secondary = as.numeric(out$secondary),
    time_s    = timeS,
    iters     = as.integer(out$iters)
  )
}
