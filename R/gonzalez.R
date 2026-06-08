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
    mean(sub[lower.tri(sub)])  # nocov
  }
}

#' Deterministic Gonzalez furthest-point selection
#'
#' Greedy k-centre selection (Gonzalez 1985). Iteratively selects the point
#' furthest from the current selection, a 2-approximation to the k-centre
#' problem. The quality of the result depends on the first (seed) point; by
#' default `Gonzalez()` runs an **ensemble** of cheap `O(N)` seeding strategies
#' and keeps the selection with the largest minimum pairwise distance
#' ([TkScore()]). The default ensemble is the two deterministic `O(N)` seeds
#' `"centroid"` and `"peripheral"` together with `n_random` (3) reproducible
#' `"random_furthest"` starts -- a best-of-five selection. The `"centroid"`
#' seed is computed from coordinates, so on a distance matrix the default drops
#' it and keeps `"peripheral"` plus the random starts. Costlier `O(N^2)`
#' anchors (`"diameter"`, `"anti_medoid"`, `"medoid"`, `"rowsum"`, `"rownorm"`)
#' are available as opt-in `seed` strategies.
#'
#' `Gonzalez()` accepts the distances in whichever of three forms suits the
#' data, all returning identical selections on the same metric:
#' \describe{
#'   \item{a **distance matrix** (`d`)}{a `dist` object or square matrix, held
#'     in full;}
#'   \item{a **coordinate matrix** (`points`)}{each needed distance is
#'     recomputed from coordinates on the fly in `O(N)` memory, never
#'     materialising the `N x N` matrix (Euclidean data only);}
#'   \item{a **distance-column oracle** (a function passed as `d`)}{for metrics
#'     with neither a stored matrix nor a coordinate embedding -- e.g.
#'     tree-to-tree distances computed on demand. See *Distance-column oracle*.}
#' }
#'
#' @section Distance-column oracle:
#' When `d` is a function it is treated as a closure `colFn(i)` returning, for a
#' single 1-based index `i`, the length-`N` vector of distances from element `i`
#' to every element. Gonzalez calls it once per selected element, maintaining a
#' running nearest-distance vector, so the `N x N` matrix is never built:
#' `O(N * n)` oracle calls and `O(N)` memory. The self-distance at position `i`
#' may take any non-negative value; it is masked before use. Because the count
#' of elements cannot be inferred from the closure, `N` must be supplied.
#' Only an integer `seed` (a `first` index) or a single deterministic two-sweep
#' peripheral seed is reachable from an oracle; the other anchors need either
#' coordinates or `O(N^2)` work. An explicitly named or ensemble `seed` on this
#' path is ignored with a warning and the peripheral seed is used instead.
#'
#' @param d A `dist` object, a square symmetric numeric matrix of pairwise
#'   distances, or a **distance-column oracle** function (see
#'   *Distance-column oracle*). Ignored when `points` is supplied.
#' @param n Integer: number of points to select. If `n >= N`, all
#'   indices are returned.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the selection is computed directly from coordinates in `O(N * n * dim)`
#'   time and `O(N)` memory, never materialising the `N x N` distance matrix
#'   (`d` is then unused). For Euclidean data the returned indices are
#'   identical to the matrix path. Only complete (non-`NA`) data is supported.
#' @param N Integer: the total number of elements. Required (and used) only on
#'   the distance-column oracle path, where it cannot be inferred from the
#'   closure; ignored for the matrix and coordinate paths.
#' @param progress Logical; show a progress bar during greedy selection on the
#'   distance-column oracle path (the only path slow enough to warrant one).
#'   Default: `TRUE` in interactive sessions, `FALSE` otherwise
#'   (`getOption("MaxMin.progress", interactive())`).
#' @param seed Integer or character (scalar or vector). An **integer** gives the
#'   explicit 1-based index of the first selected point (a single bare Gonzalez
#'   pass). A **length-1 character** names a single seeding strategy:
#'   `"centroid"` (coordinates only), `"peripheral"` (two-sweep
#'   diameter-endpoint approximation), `"random_furthest"` (furthest point from
#'   a fixed-seed random pivot), `"diameter"`, `"anti_medoid"`, `"medoid"`,
#'   `"rowsum"`, `"rownorm"`, or `"first"` (index 1). A **length > 1 character
#'   vector** requests an ensemble: each named anchor runs a full Gonzalez pass
#'   and the best result by [TkScore()] is returned with `strategy_results` and
#'   `winning_strategy` (character vector of all tied-best strategies)
#'   attributes. The `"random_furthest"` token expands to `n_random` starts,
#'   labelled `random_furthest1`, `random_furthest2`, ... Valid ensemble anchors:
#'   any subset of `c("centroid", "peripheral", "random_furthest", "diameter",
#'   "anti_medoid", "medoid", "rowsum", "rownorm")` (`"centroid"` requires
#'   `points`). Default: `c("centroid", "peripheral", "random_furthest")`. See
#'   [MaxMinSeed()] for anchor definitions. On the distance-column oracle path
#'   only an integer `seed` is honoured; a named or ensemble `seed` there warns
#'   and falls back to the peripheral seed (see *Distance-column oracle*).
#' @param n_random Integer; the number of fixed-seed random-furthest starts the
#'   `"random_furthest"` ensemble token expands to. The random pivots are drawn
#'   from a fixed internal seed, isolated from the ambient RNG, so the default
#'   selection is reproducible across sessions and machines. Default `3`; `0`
#'   contributes no random starts.
#' @return Integer vector of length `min(n, N)` of selected indices.
#' @seealso [MaxMinSeed()] for the seed indices alone; [DropAddTS()] and
#'   [ExactMaxMin()] for higher-effort solvers.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' # Default: best of the O(N) seeds (centroid, peripheral, 3 random-furthest):
#' Gonzalez(d, 5L)
#' # Raise the number of random-furthest starts:
#' Gonzalez(d, 5L, n_random = 8L)
#' # Custom two-anchor ensemble:
#' Gonzalez(d, 5L, seed = c("diameter", "anti_medoid"))
#' # A single strategy:
#' Gonzalez(d, 5L, seed = "diameter")
#' # An explicit start index (integer seed):
#' Gonzalez(d, 5L, seed = 1L)
#' # Matrix-free coordinate path (identical result, O(N) memory):
#' Gonzalez(n = 5L, points = pts, seed = 1L)
#' # Distance-column oracle: supply one column at a time, never the full matrix.
#' dm <- as.matrix(d)
#' colFn <- function(i) dm[, i]
#' identical(Gonzalez(colFn, 5L, N = nrow(dm), seed = 1L),
#'           Gonzalez(d, 5L, seed = 1L))
#' @export
Gonzalez <- function(d = NULL, n,
                     seed = .kDefaultEnsemble, n_random = 3L,
                     points = NULL, N = NULL,
                     progress = getOption("MaxMin.progress", interactive())) {
  seedMissing <- missing(seed)
  n_random <- as.integer(n_random)
  if (length(n_random) != 1L || is.na(n_random) || n_random < 0L) {
    stop("`n_random` must be a single non-negative integer")
  }
  if (is.numeric(seed) || is.integer(seed)) {
    first <- as.integer(seed)
    seed  <- "first"
  } else if (length(seed) > 1L) {
    first <- NULL
  } else {
    first <- NULL
    seed  <- match.arg(seed, choices = c("centroid", "peripheral",
                                         "random_furthest", "diameter",
                                         "anti_medoid", "medoid", "rowsum",
                                         "rownorm", "first"))
  }

  # Distance-column oracle path: `d` is a closure returning one matrix column
  # at a time, for metrics with neither a stored matrix nor a coordinate
  # embedding (e.g. on-demand tree-to-tree distances). The selection is
  # identical to the matrix path given the same `first`; only an integer
  # `seed` (a `first` index) or the deterministic peripheral seed is reachable
  # here, since the richer anchors need O(N^2) work (see Details).
  if (is.function(d)) {
    # A named/character `seed` is unreachable from an oracle (it would need the
    # whole matrix); warn rather than silently substituting the peripheral seed.
    # `first` is non-NULL only for an integer `seed`, which *is* honoured.
    if (!seedMissing && is.null(first)) {
      warning("distance-column oracle path: only an integer `seed` (a `first` ",
              "index) is honoured; using the deterministic peripheral seed")
    }
    return(.GonzalezColumn(colFn = d, N = N, n = n, first = first,
                           progress = progress))
  }

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

  if (length(seed) > 1L) {
    if (usePoints) {
      return(.GonzEnsembleFromPoints(points, n, seed, n_random))
    }
    # Matrix path: `"centroid"` is coordinate-only, so drop it (the remaining
    # O(N) seeds cover that role here). Warn only if it was named explicitly,
    # not when filtering the default ensemble.
    matrix_anchors <- seed[seed != "centroid"]
    if (!seedMissing && length(matrix_anchors) < length(seed)) {
      warning("`centroid` seed requires coordinates; it is dropped on the ",
              "distance-matrix path, where `peripheral` covers the same role")
    }
    return(.GonzEnsemble(d, n, matrix_anchors, n_random))
  }

  if (!usePoints && seed == "centroid") {
    stop("`centroid` seed requires coordinates; supply `points=` or use ",
         "`peripheral` on the distance-matrix path")
  }
  s <- if (usePoints) .MaxMinSeedPoints(points, seed) else .MaxMinSeed(d, seed)
  greedy(s)
}

#' Gonzalez maximin from a distance-column oracle (worker)
#'
#' Implements the distance-column oracle path of [Gonzalez()] (dispatched there
#' when `d` is a function); see that function's *Distance-column oracle* section
#' for the user-facing contract. At each greedy step the distances from the
#' newly selected element to all `N` elements are obtained from `colFn`, and a
#' running nearest-distance vector is maintained, so the `N x N` distance matrix
#' is never materialised: `O(N * n)` oracle calls and `O(N)` memory.
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
#' @return Integer vector of length `min(n, N)` of selected indices.
#' @keywords internal
.GonzalezColumn <- function(colFn, N, n, first = NULL,
                            progress = getOption("MaxMin.progress", interactive())) {
  if (!is.function(colFn)) {
    stop("`colFn` must be a function of one index returning numeric(N)")  # nocov
  }
  if (is.null(N)) {
    stop("`N` (the element count) must be supplied when `d` is a ",
         "distance-column function")
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
#' @param colFn Column oracle; see [Gonzalez()].
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
    .pb <- cli::cli_progress_bar("Gonzalez (column oracle)", total = n - 1L)
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
#' @param colFn Column oracle; see [Gonzalez()].
#' @param N Integer element count.
#' @return Integer index of the seed.
#' @keywords internal
.PeripheralSeedColumn <- function(colFn, N) {
  s1 <- which.max(as.numeric(colFn(1L)))
  which.max(as.numeric(colFn(s1)))
}
