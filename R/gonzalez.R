# gonzalez.R

#' Coerce distance input to a square matrix, skipping the round-trip when
#' already a matrix.
#' @param d A `dist` object or a square symmetric numeric matrix.
#' @return A square numeric matrix.
#' @keywords internal
.AsDistMatrix <- function(d) {
  if (inherits(d, "dist")) return(as.matrix(d))
  if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a `dist` object or a square numeric matrix")
  }
  d
}

#' Gonzalez maximin from a single starting index
#'
#' Internal helper: greedy furthest-point selection starting from a specified
#' index.
#'
#' @param d Square pairwise distance matrix.
#' @param n Integer: target subsample size (`>= 1`).
#' @param first Integer: index of the first selected point.
#' @return Integer vector of length `n` of selected row/col indices.
#' @keywords internal
.MaximinFrom <- function(d, n, first) {
  MaximinFrom_cpp(d, as.integer(n), as.integer(first))
}

#' Coerce coordinate input for the on-the-fly (matrix-free) samplers
#'
#' The coordinate paths require a complete numeric `N x dim` matrix with
#' `double` storage; the C++ kernels reproduce `stats::dist()`'s exact
#' Euclidean bits, which is only defined for complete data.
#' @param points A numeric matrix (or coercible) of point coordinates.
#' @return A `double` numeric matrix.
#' @keywords internal
.AsPointsMatrix <- function(points) {
  if (!is.matrix(points)) {
    points <- as.matrix(points)
  }
  if (!is.numeric(points)) {
    stop("`points` must be a numeric coordinate matrix")
  }
  if (storage.mode(points) != "double") {
    storage.mode(points) <- "double"
  }
  if (anyNA(points)) {
    stop("`points` must not contain NA; the coordinate path is for complete ",
         "Euclidean data")
  }
  points
}

#' Gonzalez maximin from coordinates (matrix-free)
#'
#' Coordinate counterpart of [.MaximinFrom()]: greedy furthest-point selection
#' that recomputes each needed distance column from `points` on the fly,
#' never materialising the `N x N` matrix. Bit-identical selection to the
#' matrix path on Euclidean data.
#'
#' @param points A `double` `N x dim` coordinate matrix.
#' @param n Integer subsample size.
#' @param first Integer index of the first selected point.
#' @param mask Integer 1-based index of a point to forbid from selection
#'   (`0L` = none); used by the anti-medoid path to exclude the medoid.
#' @return Integer vector of selected indices.
#' @keywords internal
.MaximinFromPoints <- function(points, n, first, mask = 0L) {
  MaximinFromPoints_cpp(points, as.integer(n), as.integer(first),
                        as.integer(mask))
}

#' Minimum pairwise distance within a selection, from coordinates
#'
#' Coordinate counterpart of `.SubsetScore(d, idx, "min_pairwise")`. Computes
#' `stats::dist()` on the selected sub-coordinates only (`k x k`, never the
#' full matrix); the per-pair bits are identical to the corresponding entries
#' of the full distance matrix, so the returned scalar matches the matrix path.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param idx Integer indices of the selection.
#' @return Numeric scalar; `NA_real_` if `length(idx) < 2`.
#' @keywords internal
.MinPairwiseFromPoints <- function(points, idx) {
  if (length(idx) < 2L) {
    return(NA_real_)
  }
  min(stats::dist(points[idx, , drop = FALSE]))
}

#' Score a Gonzalez subset by its minimum (or mean) pairwise distance
#'
#' @param d Full pairwise distance matrix.
#' @param idx Integer indices of selected rows/cols.
#' @param objective `"min_pairwise"` (default; canonical Gonzalez T_k) or
#'   `"mean_pairwise"`.
#' @return Numeric scalar; `NA` if `length(idx) < 2`.
#' @keywords internal
.SubsetScore <- function(d, idx, objective = c("min_pairwise", "mean_pairwise")) {
  objective <- match.arg(objective)
  if (length(idx) < 2L) {
    return(NA_real_)
  }
  sub <- d[idx, idx, drop = FALSE]
  if (objective == "min_pairwise") {
    diag(sub) <- Inf
    min(sub)
  } else {
    mean(sub[lower.tri(sub)])
  }
}

#' Deterministic Gonzalez furthest-point selection
#'
#' Greedy k-centre selection (Gonzalez 1985). Iteratively selects the point
#' furthest from the current selection, a 2-approximation to the k-centre
#' problem. The quality of the result depends on the first (seed) point; by
#' default `Gonzalez()` runs an **ensemble** of cheap deterministic peripheral
#' seeding strategies and keeps the selection with the largest minimum pairwise
#' distance ([TkScore()]). A single strategy may be requested through `seed`,
#' or an explicit start index through `first`.
#'
#' @param d A `dist` object or a square symmetric numeric matrix of pairwise
#'   distances. Ignored when `points` is supplied.
#' @param n Integer: number of points to select. If `n >= nrow(d)`, all
#'   indices are returned.
#' @param first Integer: explicit index of the first selected point. When
#'   supplied (non-`NULL`) it overrides `seed`, giving a single bare Gonzalez
#'   pass from that index. Default `NULL`.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the selection is computed directly from coordinates in `O(N * n * dim)`
#'   time and `O(N)` memory, never materialising the `N x N` distance matrix
#'   (`d` is then unused). For Euclidean data the returned indices are
#'   identical to the matrix path. Only complete (non-`NA`) data is supported.
#' @param seed Character: the seeding strategy used when `first` is `NULL`. One
#'   of `"ensemble"` (default; run the four deterministic anchors below and keep
#'   the best by [TkScore()]), `"diameter"`, `"anti_medoid"`, `"medoid"`,
#'   `"rowsum"`, `"rownorm"`, `"peripheral"` (a two-sweep diameter-endpoint
#'   approximation), or `"first"` (index 1, a bare Gonzalez pass). See
#'   [MaxMinSeed()] for the anchor definitions.
#' @param anchors Character vector of anchor names used only when
#'   `seed = "ensemble"`: any subset of `c("diameter", "anti_medoid", "rowsum",
#'   "rownorm")`. Default: all four. The returned vector then carries
#'   `strategy_results` and `winning_strategy` attributes.
#' @return Integer vector of length `min(n, nrow(d))` of selected indices.
#' @seealso [MaxMinSeed()] for the seed indices alone; [GonzalezColumn()] for a
#'   distance-column oracle; [DropAddTS()] and [ExactMaxMin()] for higher-effort
#'   solvers.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' # Ensemble of peripheral seeds (default):
#' Gonzalez(d, 5L)
#' # A single strategy:
#' Gonzalez(d, 5L, seed = "diameter")
#' # An explicit start index:
#' Gonzalez(d, 5L, first = 1L)
#' # Matrix-free coordinate path (identical result, O(N) memory):
#' Gonzalez(n = 5L, points = pts, first = 1L)
#' @export
Gonzalez <- function(d = NULL, n, first = NULL, points = NULL,
                     seed = c("ensemble", "diameter", "anti_medoid", "medoid",
                              "rowsum", "rownorm", "peripheral", "first"),
                     anchors = c("diameter", "anti_medoid", "rowsum",
                                 "rownorm")) {
  seed <- match.arg(seed)
  usePoints <- !is.null(points)
  if (usePoints) {
    points <- .AsPointsMatrix(points)
    nPts <- nrow(points)
  } else {
    d <- .AsDistMatrix(d)
    nPts <- nrow(d)
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
  }
  if (n >= nPts) return(seq_len(nPts))
  if (n == 0L)   return(integer(0))

  greedy <- if (usePoints) {
    function(s) .MaximinFromPoints(points, n, first = as.integer(s))
  } else {
    function(s) .MaximinFrom(d, n, first = as.integer(s))
  }

  # An explicit `first` is a single bare Gonzalez pass, overriding `seed`.
  if (!is.null(first)) {
    return(greedy(first))
  }

  if (seed == "ensemble") {
    return(if (usePoints) {
      .GonzEnsembleFromPoints(points, n, anchors)
    } else {
      .GonzEnsemble(d, n, anchors)
    })
  }

  s <- if (usePoints) .MaxMinSeedPoints(points, seed) else .MaxMinSeed(d, seed)
  greedy(s)
}

#' Gonzalez maximin from a distance-column oracle (matrix-free, any metric)
#'
#' Closure-driven counterpart of [Gonzalez()] for spaces with no coordinate
#' embedding (e.g. phylogenetic trees), where neither the matrix path nor the
#' Euclidean coordinate path applies. At each greedy step the distances from
#' the newly selected element to all `N` elements are obtained from `colFn`,
#' and a running nearest-distance vector is maintained. The `N x N` distance
#' matrix is never materialised: `O(N * n)` oracle calls and `O(N)` memory.
#'
#' On a symmetric metric the selection is identical to
#' `Gonzalez(d, n, first = first)` given the same starting index, because
#' `colFn(i)` returns column `i` of the distance matrix.
#'
#' Only the deterministic two-sweep peripheral seed is available here: the
#' richer anchors of [Gonzalez()] (diameter, anti-medoid, row-sum, row-norm)
#' need `O(N^2)` work and are unreachable from an `O(N)`-per-call oracle.
#'
#' @param colFn A function of a single 1-based index `i` returning a length-`N`
#'   numeric vector of distances from element `i` to every element. The
#'   self-distance at position `i` may take any non-negative value; it is
#'   masked before use.
#' @param N Integer: the total number of elements. It cannot be inferred from
#'   `colFn`, so it must be supplied.
#' @param n Integer: number of elements to select. If `n >= N`, all indices are
#'   returned.
#' @param first Integer index of the first selected element, or `NULL`
#'   (default) to use a deterministic peripheral seed computed from two oracle
#'   sweeps: the element furthest from element 1, then the element furthest
#'   from that (a diameter-endpoint approximation).
#' @param progress Logical; show a progress bar during greedy selection.
#'   Default: `TRUE` in interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
#' @return Integer vector of length `min(n, N)` of selected indices.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- as.matrix(dist(pts))
#' colFn <- function(i) d[, i]
#' identical(GonzalezColumn(colFn, nrow(d), 5L, first = 1L),
#'           Gonzalez(d, 5L, first = 1L))
#' @seealso [Gonzalez()] for the matrix and coordinate paths.
#' @export
GonzalezColumn <- function(colFn, N, n, first = NULL,
                            progress = getOption("MaxMin.progress", interactive())) {
  if (!is.function(colFn)) {
    stop("`colFn` must be a function of one index returning numeric(N)")
  }
  N <- as.integer(N)
  n <- as.integer(n)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    stop("`N` must be a single positive integer")
  }
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
  }
  if (n >= N) return(seq_len(N))
  if (n == 0L) return(integer(0))
  if (is.null(first)) {
    first <- .PeripheralSeedColumn(colFn, N)
  }
  first <- as.integer(first)
  if (length(first) != 1L || is.na(first) || first < 1L || first > N) {
    stop("`first` must be a single index in [1, N]")
  }
  if (n == 1L) return(first)
  .MaximinFromColumn(colFn, N, n, first, progress = progress)
}

#' Gonzalez maximin from a distance-column oracle (worker)
#'
#' Mirrors `MaximinFrom_cpp()` (src/maximin.cpp), substituting an on-demand
#' `colFn(i)` call for the matrix-column read `d[, i]`. `which.max()` uses R's
#' first-maximum (strict `>`) rule, matching the kernel's tie-breaking, so the
#' selection is identical to the matrix path on symmetric input.
#' @param colFn Column oracle; see [GonzalezColumn()].
#' @param N Integer element count.
#' @param n Integer subset size (`>= 2`).
#' @param first Integer seed index.
#' @return Integer vector of selected indices.
#' @keywords internal
.MaximinFromColumn <- function(colFn, N, n, first, progress = FALSE) {
  selected <- integer(n)
  selected[1L] <- first
  min_dist <- as.numeric(colFn(first))
  if (length(min_dist) != N) {
    stop("`colFn` must return a numeric vector of length N = ", N,
         "; got length ", length(min_dist))
  }
  min_dist[first] <- -Inf                 # mask seed before the loop
  if (progress) {
    .pb <- cli::cli_progress_bar("GonzalezColumn", total = n - 1L)
  }
  for (k in seq_len(n - 1L) + 1L) {
    best <- which.max(min_dist)           # first global max (ties -> first)
    selected[k] <- best
    min_dist[best] <- -Inf                # mask before pmin so self-dist 0
                                          # cannot overwrite -Inf
    min_dist <- pmin.int(min_dist, as.numeric(colFn(best)))
    if (progress) cli::cli_progress_update(id = .pb)
  }
  selected
}

#' Deterministic peripheral seed from a column oracle
#'
#' Two oracle sweeps, no RNG: the element furthest from element 1, then the
#' element furthest from that. The second is a diameter-endpoint approximation
#' and a markedly better Gonzalez anchor than an arbitrary start, at the cost
#' of two of the `O(N * n)` sweeps. The richer peripheral anchors (diameter,
#' anti-medoid) need `O(N^2)` work and are unreachable from a column oracle.
#' @param colFn Column oracle; see [GonzalezColumn()].
#' @param N Integer element count.
#' @return Integer index of the seed.
#' @keywords internal
.PeripheralSeedColumn <- function(colFn, N) {
  s1 <- which.max(as.numeric(colFn(1L)))
  which.max(as.numeric(colFn(s1)))
}
