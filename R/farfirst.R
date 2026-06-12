# farfirst.R

#' Coerce distance input to a square matrix, skipping the round-trip when
#' already a matrix.
#' @param d A `dist` object or a square numeric matrix.
#' @return A square numeric matrix.
#' @details
#' Symmetry is not checked; an `O(N^2)` check is intentionally omitted.
#' Asymmetric matrices are silently accepted, and the algorithm treats
#' \eqn{d_{ij}} and \eqn{d_{ji}} as independent values.
#' @keywords internal
.AsDistMatrix <- function(d) {
  if (inherits(d, "dist")) {
    d <- as.matrix(d)
  } else if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a `dist` object or a square numeric matrix")
  }
  # NA/NaN/Inf propagate silently through pmin.int()/which.max() and can yield a
  # repeated index in the selection; reject them up front (the coordinate path
  # already errors via .AsPointsMatrix()'s anyNA check).
  if (anyNA(d) || any(!is.finite(d))) {
    stop("distance matrix must not contain NA/NaN/Inf")
  }
  d
}

#' Gonzalez maximin from a single starting index
#'
#' Internal helper: greedy furthest-point selection starting from a specified
#' index.
#'
#' @param d Square pairwise distance matrix.
#' @param k Integer: target subsample size (`>= 1`).
#' @param first Integer: index of the first selected point.
#' @return Integer vector of length `k` of selected row/col indices.
#' @keywords internal
.MaximinFrom <- function(d, k, first) {
  MaximinFrom_cpp(d, as.integer(k), as.integer(first))
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
#' @inheritParams .MaximinFrom
#' @param mask Integer 1-based index of a point to forbid from selection
#'   (`0L` = none); used by the anti-medoid path to exclude the medoid.
#' @return Integer vector of selected indices.
#' @keywords internal
.MaximinFromPoints <- function(points, k, first, mask = 0L) {
  MaximinFromPoints_cpp(points, as.integer(k), as.integer(first),
                        as.integer(mask))
}

#' Promote the kernel's free `t_k` score to the user-facing `score` attribute
#'
#' The maximin kernels attach the selection's minimum pairwise distance as a
#' `t_k` attribute (computed during the greedy pass at no extra cost). The
#' ensemble drivers read it via [base::attr()] into `strategy_results`; a bare
#' single pass exposes it directly as the `score` attribute, matching
#' [DropAdd()] and [Grasp()].
#' @param idx Integer vector returned by a maximin kernel.
#' @return `idx` with its `t_k` attribute renamed to `score`.
#' @keywords internal
.PromoteScore <- function(idx) {
  attr(idx, "score") <- attr(idx, "t_k")
  attr(idx, "t_k")   <- NULL
  idx
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
#' Greedy _k_-centre selection \insertCite{Gonzalez1985,Hochbaum1985}{MaxMin}.
#' Iteratively selects the point furthest from the current selection, a
#' 2-approximation to the _k_-centre problem.
#' The quality of the result depends on the first (seed) point; by default
#' `FarFirst()` runs three starts from randomly selected peripheral seeds.
#' The deterministic \eqn{O(N)} anchors (`"centroid"`, `"peripheral"`) and the
#' costlier \eqn{O(N^2)} anchors (`"diameter"`, `"anti_medoid"`,
#' `"rowsum"`, `"rownorm"`) are alternative `seed` strategies.
#'
#' Distances may be provided as:
#' \describe{
#'   \item{a **distance matrix** (`d`)}{a `dist` object or square matrix, held
#'     in full;}
#'   \item{a **coordinate matrix** (`points`)}{each needed distance is
#'     recomputed from coordinates on the fly in `O(N)` memory
#'     (Euclidean data only);}
#'   \item{a distance function (function passed as `d`)}{for metrics
#'     with neither a stored matrix nor a coordinate embedding, where distances
#'     may be computed on demand.}
#' }
#'
#' @section Distance function:
#' When `d` is a function, it will be passed a single 1-based index `i`, and
#' should return the distances from element `i` to every element in turn,
#' optionally omitting entry `i`, the self-distance.
#' The function will be called once per selected element, to avoid building a
#' complete `N x N` matrix.
#'
#' @param d A `dist` object, a square numeric matrix of pairwise distances, or
#'   a distance function (see *§Distance function*). Asymmetric matrices are
#'   accepted; the algorithm treats \eqn{d_{ij}} and \eqn{d_{ji}} as
#'   independent. Ignored when `points` is supplied.
#' @param k Integer: number of points to select.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the selection is computed directly from coordinates in `O(N * k * dim)`
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
#' @param method Integer or character (scalar or vector); how to seed the
#'   greedy pass (matching the `method` argument of [MaxMinSeed()]). An
#'   **integer** gives the explicit 1-based index of the first selected point
#'   (a single bare Gonzalez pass). A **length-1 character** names a single
#'   deterministic seeding strategy run as one bare pass: `"centroid"`
#'   (coordinates only), `"peripheral"` (two-sweep diameter-endpoint
#'   approximation), `"diameter"`, `"anti_medoid"`, `"rowsum"`,
#'   `"rownorm"`, or `"first"` (index 1). A **length > 1 character vector** --
#'   or the lone `"random_furthest"` token -- requests an ensemble: each named
#'   anchor runs a full Gonzalez pass and the best result by [MinDist()] is
#'   returned with `strategy_results` and `winning_strategy` (character vector
#'   of all tied-best strategies) attributes. The `"random_furthest"` token
#'   expands to one start per element of `pivots`, labelled `random_furthest1`,
#'   `random_furthest2`, ...; named on its own it still runs the ensemble (one
#'   pass per pivot), so a single random start is best obtained via
#'   [MaxMinSeed()]. Valid ensemble anchors: any subset of `c("centroid",
#'   "peripheral", "random_furthest", "diameter", "anti_medoid", "rowsum",
#'    "rownorm")` (`"centroid"` requires `points`). Default:
#'   `"random_furthest"` (three random starts; see `pivots`). See [MaxMinSeed()]
#'   for anchor definitions. On the distance-column oracle path only an integer
#'   `method` is honoured; a named or ensemble `method` there warns and falls
#'   back to the peripheral seed (see *Distance-column oracle*).
#' @param pivots Integer vector of pivot indices over which the
#'   `"random_furthest"` ensemble token expands: each pivot contributes one
#'   start, seeded at the point furthest from it, so the vector's length sets
#'   the number of random-furthest starts. Left unspecified, three pivots are
#'   drawn with the session RNG (`sample.int(N, 3)`; set a seed for a
#'   reproducible selection). Pass `integer(0)`, `NA`, or `NULL` to disable the
#'   random starts, or an index vector to choose the pivots (and their count)
#'   explicitly. Disabling the random starts errors under the default `method`
#'   (which names only `"random_furthest"`, leaving no anchor); pair it with a
#'   deterministic `method` such as `"peripheral"`.
#' @param nseeds Optional integer: run a distinct-seed random restart. Random
#'   pivots are drawn with the session RNG and each one's furthest-point seed is
#'   collected, de-duplicated, until `nseeds` *distinct* seeds are found (or the
#'   reachable pool is exhausted); Gonzalez runs from each and the best \eqn{T_k}
#'   is returned, with `strategy_results` / `winning_strategy` labelled
#'   `random_furthest1`, `random_furthest2`, ... This is the "give a count, not a
#'   list" counterpart to `pivots`: where `pivots` runs one start per supplied
#'   index, `nseeds` searches for that many *distinct* peripheral seeds, never
#'   wasting a Gonzalez pass on a duplicate. Set a seed (`set.seed()`) for a
#'   reproducible selection. When supplied, `nseeds` overrides `method` and
#'   `pivots` (a warning is issued if either was also set explicitly); it is not
#'   available on the distance-column oracle path. Default `NULL` (use `method`).
#' @return Integer vector of length `min(k, N)` of selected indices, in
#'   the order they were selected.
#'   The achieved \eqn{T_k} (the selection's minimum pairwise
#'   distance) is attached as attribute `score`. An ensemble `method`
#'   additionally carries `strategy_results` and `winning_strategy` attributes.
#'   The vector has class `"MaxMinSelection"` and prints as a one-line summary
#'   (see [print.MaxMinSelection]); it is otherwise an ordinary integer vector
#'   and indexes a matrix or coordinate set directly.
#' @references \insertAllCited{}
#' @seealso [MaxMinSeed()] for the seed indices alone; [DropAdd()] and
#'   [ExactMaxMin()] for higher-effort solvers.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' # Default: best of three random-furthest starts (set.seed for reproducibility):
#' FarFirst(d, 5L)
#' # More random-furthest starts (length of `pivots` sets the count):
#' FarFirst(d, 5L, pivots = sample.int(nrow(as.matrix(d)), 8L))
#' # Or choose the pivots explicitly:
#' FarFirst(d, 5L, pivots = c(1L, 10L, 20L))
#' # Custom two-anchor ensemble:
#' FarFirst(d, 5L, method = c("diameter", "anti_medoid"))
#' # A single strategy:
#' FarFirst(d, 5L, method = "diameter")
#' # An explicit start index (integer method):
#' FarFirst(d, 5L, method = 1L)
#' # Matrix-free coordinate path (identical result, O(N) memory):
#' FarFirst(k = 5L, points = pts, method = 1L)
#'
#' # Distance-column oracle: supply one column at a time, never the full matrix.
#' data("USArrests")
#' arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
#' StateDist <- function(i) {
#'   diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
#'   sqrt(rowSums(diffs ^ 2))
#' }
#' idx <- FarFirst(StateDist, k = 4L, N = nrow(arrestTypes), method = 1L)
#' arrestTypes[idx, ]
#' @export
FarFirst <- function(d = NULL, k,
                     method = .kDefaultEnsemble, pivots = NULL, nseeds = NULL,
                     points = NULL, N = NULL,
                     progress = getOption("MaxMin.progress", interactive())) {
  methodMissing <- missing(method)
  pivotsMissing <- missing(pivots)

  # Every non-empty FarFirst result is a `MaxMinSelection` (a self-describing
  # integer index vector); see print.MaxMinSelection(). Stamping it at each
  # exit keeps the matrix, coordinate and oracle paths byte-identical.
  Classify <- function(x) .AsMaxMinSelection(x, "FarFirst")

  # All single-strategy names (the ensemble anchors plus the single-pass-only
  # `"first"`); .kPointEnsembleSeeds is the ensemble-eligible subset.
  validMethods <- c(.kPointEnsembleSeeds, "first")

  if (is.numeric(method)) {
    # An integer `method` is the explicit 1-based first index (a single bare
    # Gonzalez pass). Guard it here: NA/NaN would otherwise reach the C++
    # kernel as INT_MIN, and a zero-length or multi-element value as an opaque
    # Rcpp "Expecting a single value" error.
    if (length(method) != 1L || !is.finite(method)) {
      stop("an integer `method` (the first index) must be a single finite value")
    }
    first  <- as.integer(method)
    method <- "first"
  } else if (length(method) > 1L) {
    # A multi-element `method` requests an ensemble. The downstream
    # match.arg(several.ok = TRUE) silently *drops* names that fail to match,
    # so validate explicitly here against the ensemble anchors (centroid
    # included -- the matrix path drops it later with a warning, not an error).
    bad <- setdiff(method, .kPointEnsembleSeeds)
    if (length(bad)) {
      stop("unknown seed/method strateg", if (length(bad) > 1L) "ies: " else "y: ",
           paste(bad, collapse = ", "))
    }
    first  <- NULL
  } else {
    first  <- NULL
    method <- match.arg(method, choices = validMethods)
  }

  # Distance-column oracle path: `d` is a closure returning one matrix column
  # at a time, for metrics with neither a stored matrix nor a coordinate
  # embedding (e.g. on-demand tree-to-tree distances). The selection is
  # identical to the matrix path given the same `first`; only an integer
  # `method` (a `first` index) or the deterministic peripheral seed is reachable
  # here, since the richer anchors need O(N^2) work (see Details).
  if (is.function(d)) {
    if (!is.null(nseeds)) {
      stop("`nseeds` (distinct-seed random restart) is not supported on the ",
           "distance-column oracle path; supply an integer `method` (a `first` ",
           "index) instead")
    }
    # A named/character `method` is unreachable from an oracle (it would need
    # the whole matrix); warn rather than silently substituting the peripheral
    # seed. `first` is non-NULL only for an integer `method`, which *is* honoured.
    if (!methodMissing && is.null(first)) {
      warning("distance-column oracle path: only an integer `method` (a `first` ",
              "index) is honoured; using the deterministic peripheral seed")
    }
    return(Classify(.GonzalezColumn(colFn = d, N = N, k = k, first = first,
                                    progress = progress)))
  }

  usePoints <- !is.null(points)
  if (usePoints) {
    points <- .AsPointsMatrix(points)
    nPts <- nrow(points)
  } else {
    d <- .AsDistMatrix(d)
    nPts <- nrow(d)
  }
  # Pre-check before as.integer(): a finite, length-1, non-negative value.
  # Pre-checking avoids the spurious base-R coercion warning on `k = Inf`.
  if (length(k) != 1L || !is.finite(k) || k < 0L) {
    stop("`k` must be a single non-negative integer")
  }
  k <- min(as.integer(k), nPts)
  if (k == 0L) return(integer(0))

  # Distinct-seed random restart: draw random pivots, take each one's
  # furthest-point seed, de-duplicate, until `nseeds` distinct seeds are
  # collected (or the reachable pool is exhausted); run Gonzalez from each and
  # keep the best. The "give a count" counterpart to an explicit `pivots` list;
  # it overrides `method`/`pivots`. RNG is the session stream (set.seed()).
  if (!is.null(nseeds)) {
    nseeds <- as.integer(nseeds)
    if (length(nseeds) != 1L || is.na(nseeds) || nseeds < 1L) {
      stop("`nseeds` must be a single positive integer")
    }
    if (!methodMissing || !pivotsMissing) {
      warning("`nseeds` runs the distinct-seed random restart and overrides ",
              "`method`/`pivots`")
    }
    seedFn <- if (usePoints) {
      function(r) which.max(EuclidColFromPoints_cpp(points, r))
    } else {
      function(r) which.max(d[, r])
    }
    seeds <- .DrawDistinctSeeds(seedFn, nPts, nseeds)
    return(Classify(if (usePoints) {
      .GonzEnsembleFromPoints(points, k, "random_furthest", pivots = seeds,
                              rfSeedFn = identity)
    } else {
      .GonzEnsemble(d, k, "random_furthest", pivots = seeds,
                    rfSeedFn = identity)
    }))
  }

  # The kernels attach a `t_k` attribute (the selection's min pairwise distance,
  # computed for free). A bare single pass exposes it as the `score` attribute
  # (matching DropAdd/Grasp); the ensemble path reads it into strategy_results.
  Greedy <- if (usePoints) {
    function(s) .PromoteScore(.MaximinFromPoints(points, k, first = as.integer(s)))
  } else {
    function(s) .PromoteScore(.MaximinFrom(d, k, first = as.integer(s)))
  }

  # An explicit `first` is a single bare Gonzalez pass, overriding `method`.
  if (!is.null(first)) {
    return(Classify(Greedy(first)))
  }

  # The ensemble path runs each named anchor as a full Gonzalez pass and keeps
  # the best by MinDist(). It is taken for a multi-anchor `method`, and also for
  # a lone `"random_furthest"` (the default): that token is inherently
  # multi-start -- it expands to one pass per `pivots` element -- so it belongs
  # here, not on the single-strategy path. (`MaxMinSeed(method =
  # "random_furthest")` still returns exactly one seed for callers who want a
  # single random start.)
  if (length(method) > 1L || "random_furthest" %in% method) {
    # Pivots for the `"random_furthest"` token. Unspecified: draw three pivots
    # with the session RNG (`set.seed()` for a reproducible selection). An empty
    # / `NA` / `NULL` `pivots` disables the random starts; a supplied index
    # vector is taken verbatim, its length setting the number of starts.
    if (pivotsMissing) {
      pivots <- if ("random_furthest" %in% method) {
        sample.int(nPts, min(.kDefaultRandomStarts, nPts))
      } else {
        integer(0)
      }
    } else if (length(pivots) == 0L || all(is.na(pivots))) {
      pivots <- integer(0)
    } else {
      pivots <- as.integer(pivots)
      if (anyNA(pivots) || any(pivots < 1L | pivots > nPts)) {
        stop("`pivots` must be indices in [1, N]")
      }
    }
    # Matrix path: `"centroid"` is coordinate-only, so drop it (the remaining
    # O(N) seeds cover that role here). Warn only if it was named explicitly,
    # not when filtering the default ensemble.
    anchors <- if (usePoints) method else method[method != "centroid"]
    if (!usePoints && !methodMissing && length(anchors) < length(method)) {
      warning("`centroid` seed requires coordinates; it is dropped on the ",
              "distance-matrix path, where `peripheral` covers the same role")
    }
    # With `"random_furthest"` the only anchor and no pivots, nothing would run.
    # This is reachable from the default `method` once the random starts are
    # disabled (`pivots = integer(0)` / `NA` / `NULL`), so fail clearly rather
    # than tripping the internal "no strategies" guard.
    if (length(setdiff(anchors, "random_furthest")) == 0L &&
        length(pivots) == 0L) {
      stop("no seed strategies to run: disabling the random-furthest starts ",
           "leaves the default ensemble with no anchor. Name a deterministic ",
           "`method` (e.g. \"peripheral\") or supply non-empty `pivots`.")
    }
    if (usePoints) {
      return(Classify(.GonzEnsembleFromPoints(points, k, anchors, pivots)))
    }
    return(Classify(.GonzEnsemble(d, k, anchors, pivots)))
  }

  if (!usePoints && method == "centroid") {
    stop("`centroid` seed requires coordinates; supply `points=` or use ",
         "`peripheral` on the distance-matrix path")
  }
  s <- if (usePoints) .MaxMinSeedPoints(points, method) else .MaxMinSeed(d, method)
  Classify(Greedy(s))
}

#' Normalise a distance-column oracle result to a masked length-`N` vector
#'
#' The user's `colFn(i)` may either report the self-distance (a length-`N`
#' vector, position `i` ignored) or omit it (a length-`N - 1` vector of the
#' distances to the other elements, in index order). Either way this returns a
#' length-`N` numeric vector with position `i` set to `-Inf`, so the downstream
#' `which.max()` / `pmin.int()` never re-select `i`. The mask invariant for the
#' oracle path lives here, not in the callers.
#' @param colFn Column oracle; see [FarFirst()].
#' @param i Integer 1-based index whose distance column is requested.
#' @param N Integer element count.
#' @return Numeric vector of length `N`, masked to `-Inf` at position `i`.
#' @keywords internal
.DistColumn <- function(colFn, i, N) {
  col <- as.numeric(colFn(i))
  len <- length(col)
  if (len == N) {
    col[i] <- -Inf                       # self reported; ignore and mask it
  } else if (len == N - 1L) {
    col <- append(col, -Inf, after = i - 1L)  # self omitted; splice the mask in
  } else {
    stop("`colFn(", i, ")` must return a numeric vector of length N = ", N,
         " (self-distance included) or N - 1 (self-distance omitted); got ",
         "length ", len)
  }
  # The self position is masked to -Inf above; any remaining NA/NaN is a genuine
  # non-self distance and would silently corrupt the greedy pass (cf. the matrix
  # path's guard in .AsDistMatrix). Reject it.
  if (anyNA(col)) {
    stop("`colFn(", i, ")` returned NA/NaN at a non-self position; ",
         "distance-column oracles must return finite distances")
  }
  col
}

#' Gonzalez maximin from a distance-column oracle
#'
#' Implements the distance-column oracle path of [FarFirst()] (dispatched there
#' when `d` is a function); see that function's *Distance-column oracle* section
#' for the user-facing contract. At each greedy step the distances from the
#' newly selected element to all `N` elements are obtained from `colFn`, and a
#' running nearest-distance vector is maintained, so the `N x N` distance matrix
#' is never materialised: `O(N * k)` oracle calls and `O(N)` memory.
#'
#' @param colFn A function that, when passed an index `i`, must return a
#' vector of distances from element `i` to either (i) every element in turn,
#' including `i`; or (ii) every other element. See [.DistColumn()].
#' @param N Integer: the total number of elements.
#' @param k Integer: number of elements to select. If `k > N`, all `N` indices
#'   are returned in Gonzalez (farthest-first) order.
#' @param first Integer index of the first selected element, or `NULL`
#'   (default) to use a deterministic peripheral seed computed from two oracle
#'   sweeps: the element furthest from element 1, then the element furthest
#'   from that (a diameter-endpoint approximation).
#' @param progress Logical; show a progress bar during greedy selection.
#' @return Integer vector of length `min(k, N)` of selected indices.
#' @keywords internal
.GonzalezColumn <- function(colFn, N, k, first = NULL,
                            progress = getOption("MaxMin.progress", interactive())) {
  if (!is.function(colFn)) {
    stop("`colFn` must be a function of one index returning numeric(N)")
  }
  if (is.null(N)) {
    stop("`N` (the element count) must be supplied when `d` is a ",
         "distance-column function")
  }
  N <- as.integer(N)
  k <- as.integer(k)
  if (length(N) != 1L || is.na(N) || N < 1L) {
    stop("`N` must be a single positive integer")
  }
  if (length(k) != 1L || is.na(k) || k < 0L) {
    stop("`k` must be a single non-negative integer")
  }
  k <- min(k, N)
  if (k == 0L) return(integer(0))
  if (is.null(first)) {
    first <- .PeripheralSeedColumn(colFn, N)
  }
  first <- as.integer(first)
  if (length(first) != 1L || is.na(first) || first < 1L || first > N) {
    stop("`first` must be a single index in [1, N]")
  }
  # A single point has no pairwise distance; report score = NA, matching the
  # matrix/coordinate kernels' k < 2 behaviour.
  if (k == 1L) return(structure(first, score = NA_real_))
  .MaximinFromColumn(colFn, N, k, first, progress = progress)
}

#' Gonzalez maximin from a distance-column oracle (worker)
#'
#' Mirrors `MaximinFrom_cpp()`, substituting an on-demand `colFn(i)` call for
#' the matrix-column read `d[, i]`. `which.max()` uses R's
#' first-maximum (strict `>`) rule, matching the kernel's tie-breaking, so the
#' selection is identical to the matrix path on symmetric input.
#' @param colFn Column oracle; see [FarFirst()].
#' @param N Integer element count.
#' @param k Integer subset size (`>= 2`).
#' @param first Integer seed index.
#' @return Integer vector of selected indices.
#' @keywords internal
.MaximinFromColumn <- function(colFn, N, k, first, progress = FALSE) {
  selected <- integer(k)
  selected[1L] <- first
  minDist <- .DistColumn(colFn, first, N)  # seed column, self already masked
  if (progress) {
    .pb <- cli::cli_progress_bar("Gonzalez (column oracle)", total = k - 1L)
  }
  # T_k = min over greedy steps of the chosen element's insertion distance
  # (minDist[best] before the pmin update), mirroring MaximinFrom_cpp's `tk` so
  # the `score` attribute matches the matrix path bit-for-bit.
  tk <- Inf
  for (i in seq_len(k - 1L) + 1L) {
    best <- which.max(minDist)           # first global max (ties -> first)
    selected[i] <- best
    if (minDist[best] < tk) tk <- minDist[best]
    # `.DistColumn` masks position `best` to -Inf, so pmin propagates the mask
    # and `best` cannot be re-selected; no explicit `minDist[best]` needed.
    minDist <- pmin.int(minDist, .DistColumn(colFn, best, N))
    if (progress) cli::cli_progress_update(id = .pb)
  }
  attr(selected, "score") <- if (is.finite(tk)) tk else NA_real_
  selected
}

#' Deterministic peripheral seed from a column oracle
#'
#' Two oracle sweeps, no RNG: the element furthest from element 1, then the
#' element furthest from that. The second is a diameter-endpoint approximation
#' and a markedly better Gonzalez anchor than an arbitrary start, at the cost
#' of two of the `O(N * n)` sweeps. The richer peripheral anchors (diameter,
#' anti-medoid) need `O(N^2)` work and are unreachable from a column oracle.
#' @param colFn Column oracle; see [FarFirst()].
#' @param N Integer element count.
#' @return Integer index of the seed.
#' @keywords internal
.PeripheralSeedColumn <- function(colFn, N) {
  s1 <- which.max(.DistColumn(colFn, 1L, N))
  which.max(.DistColumn(colFn, s1, N))
}
