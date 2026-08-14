# maxmean.R
#
# RLTS algorithm (Dieudonne et al. 2020) for the Max-Mean Dispersion Problem:
# choose a subset S (size unfixed, |S| >= 2) maximising
#   f(S) = sum_{i < j, i,j in S} d_ij / |S|
# which equals the sum of pairwise distances divided by the number of selected
# elements (not divided by the number of pairs C(|S|, 2)).
#
# The reinforcement-learning hyperparameters are fixed at the paper's tuned
# values (Section 4.1); they are not exposed, to keep the API minimal. The
# search is budgeted by wall-clock time (the paper's t_limit) and/or a hard cap
# on tabu iterations, whichever is reached first: more budget yields more
# restarts and a monotonically non-worsening best solution.

# Paper-tuned constants (Dieudonne et al. 2020, Section 4.1).
.kMaxMeanTMin    <- 5L      # minimum tabu tenure (sawtooth floor)
.kMaxMeanTMax    <- 120L    # maximum tabu tenure (T_max)
.kMaxMeanDepth   <- 50000L  # per-restart no-improve limit (alpha)
.kMaxMeanEpsilon <- 0.70    # epsilon-greedy action selection
.kMaxMeanAlpha   <- 0.50    # Q-learning rate
.kMaxMeanGamma   <- 0.50    # Q discount factor

#' RLTS tabu search for the Max-Mean Dispersion Problem
#'
#' `MaxMean()` selects a maximally dispersed subset of elements from a
#' pairwise distance matrix, maximising the *max-mean* objective:
#' \deqn{f(S) = \frac{\displaystyle\sum_{i < j,\, i,j \in S} d_{ij}}{|S|}}
#' The number of elements in the subset \eqn{|S| \ge 2} is chosen so as
#' to maximise the mean dispersion.
#'
#' `MaxMean()` implements the reinforcement-learning tabu search
#' (\acronym{RLTS}) algorithm of \insertCite{Dieudonne2020}{MaxMin}.
#' An initial solution is constructed randomly for the first restart and
#' via \eqn{Q}-learning thereafter; each initial solution is then refined by a
#' tabu search using one-flip moves (adding or removing one element per step),
#' with efficient \eqn{O(n)} neighbourhood evaluation via the contribution
#' array \eqn{P_i = \sum_{j \in S,\, j \ne i} d_{ij}}.  Restarts continue until
#' the `maxSeconds` or `maxIter` budget is reached, whichever comes first; the
#' best solution found is returned.
#'
#' The reinforcement-learning and tabu hyperparameters are fixed at the tuned
#' values reported by \insertCite{Dieudonne2020}{MaxMin} (greedy factor
#' \eqn{\epsilon = 0.7}, learning rate \eqn{\alpha = 0.5}, discount
#' \eqn{\gamma = 0.5}, maximum tabu tenure \eqn{120}, search depth 50&nbsp;000).
#'
#' @param d A `dist` object or square numeric matrix of pairwise distances;
#'   values may be negative, and an asymmetric matrix is symmetrized to
#'   \eqn{(d_{ij} + d_{ji})/2} before solving.
#' @param maxSeconds Numeric wall-clock time budget in seconds.
#' @param maxIter Numeric cap on the total tabu-search iterations across
#'   restarts; `Inf` budgets by wall-clock time alone.
#' @param useRL Logical: enable the \eqn{Q}-learning layer guiding
#'   initial-solution construction across restarts (`TRUE`, the default, costs
#'   \eqn{O(n^2)} memory; `FALSE` uses random restarts only, \eqn{O(n)} memory).
#' @templateVar progress_shows status messages are shown
#' @template progress
#'
#' @return `MaxMean()` returns an integer vector of selected 1-based indices
#'   (sorted ascending) with class `"MaxMeanSelection"` and attributes:
#'   \describe{
#'     \item{score}{numeric, achieved objective
#'       \eqn{\sum_{i<j \in S} d_{ij} / |S|}.}
#'     \item{size}{integer, number of selected elements \eqn{|S|}.}
#'     \item{time_s}{numeric, wall-clock seconds spent.}
#'     \item{iters}{numeric, total tabu-search iterations across restarts.
#'       Stored as a double, not an integer, because a long run can exceed the
#'       32-bit integer range.}
#'   }
#'   The vector has class `"MaxMeanSelection"` and prints as a one-line summary
#'   (see [print.MaxMeanSelection()]); it is otherwise an ordinary integer
#'   vector that indexes the distance matrix directly.
#'
#' @references \insertAllCited{}
#' @seealso [MeanDist()] to score an arbitrary selection under this objective;
#'   [FarFirst()], [DropAdd()] and [Grasp()] for fixed-cardinality max-min
#'   solvers.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#'
#' # Select the subset that maximises mean pairwise dispersion:
#' result <- MaxMean(d, maxSeconds = 0.001)
#' result
#' attr(result, "score")
#' @export
MaxMean <- function(d, maxSeconds = 0.1, maxIter = 1000, useRL = TRUE) {
  progress <- getOption("MaxMin.progress", interactive())

  # Asymmetry is averaged away below, so it is accepted here.
  dmat <- .AsDistMatrix(d, symmetric = FALSE)
  n    <- nrow(dmat)
  # The max-mean objective is defined for symmetric distances, but
  # .AsDistMatrix() does not enforce symmetry (the O(n^2) check is avoided as
  # elsewhere in the package). Symmetrize once here so the incremental
  # contribution array in the C++ kernel stays exact and the reported `score`
  # agrees with MeanDist(). For an already-symmetric matrix this is a
  # bit-identical no-op: (x + x) * 0.5 == x for finite x.
  dmat <- (dmat + t(dmat)) * 0.5
  # Zero the diagonal so the C++ kernel can drop its per-element `j != i` guards
  # (a no-op for the objective, which never includes a self-distance, but it
  # lets the O(n) contribution-array updates run branchless). For a `dist`
  # object or a proper distance matrix the diagonal is already zero.
  diag(dmat) <- 0

  if (!is.numeric(maxSeconds) || length(maxSeconds) != 1L ||
      is.na(maxSeconds) || maxSeconds <= 0) {
    stop("`maxSeconds` must be a single positive numeric (or Inf)")
  }
  if (!is.numeric(maxIter) || length(maxIter) != 1L ||
      is.na(maxIter) || maxIter < 1) {
    stop("`maxIter` must be a single numeric >= 1 (or Inf)")
  }
  if (!is.logical(useRL) || length(useRL) != 1L || is.na(useRL)) {
    stop("`useRL` must be a single logical (TRUE or FALSE)")
  }

  if (progress) {
    cli::cli_process_start(
      "MaxMean RLTS (n = {n}, budget = {maxSeconds}s)",
      .auto_close = FALSE
    )
  }

  t0 <- proc.time()[[3L]]

  out <- MaxMean_cpp(
    dmat          = dmat,
    time_budget_s = as.double(maxSeconds),
    iter_budget   = as.double(maxIter),
    alpha_depth   = .kMaxMeanDepth,
    T_min         = .kMaxMeanTMin,
    T_max         = .kMaxMeanTMax,
    epsilon       = .kMaxMeanEpsilon,
    alpha_rl      = .kMaxMeanAlpha,
    gamma_rl      = .kMaxMeanGamma,
    use_rl        = useRL
  )

  timeS <- proc.time()[[3L]] - t0

  if (progress) {
    kMsg    <- length(out$indices)
    fMsg    <- as.numeric(out$objective)
    iterMsg <- as.numeric(out$iters)
    cli::cli_process_done(
      msg = "MaxMean: {kMsg} elements, f = {signif(fMsg, 4)}, {iterMsg} iters, {signif(timeS, 4)}s"
    )
  }

  # `iters` is kept as a double: a long run can exceed .Machine$integer.max,
  # and as.integer() would silently return NA there.
  .AsMaxMeanSelection(structure(
    sort(as.integer(out$indices)),
    score  = as.numeric(out$objective),
    size   = length(out$indices),
    time_s = timeS,
    iters  = as.numeric(out$iters)
  ))
}
