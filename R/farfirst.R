# farfirst.R

#' Reconcile the two triangles of a distance matrix
#'
#' Internal helper: returns `d` unchanged when its triangles already agree,
#' averages them when they differ by no more than `tolerance`, and errors
#' beyond it.
#'
#' @param d A square numeric matrix, known finite.
#' @param tolerance Numeric: largest scaled discrepancy to repair.
#' @return `.Symmetrise()` returns an exactly symmetric matrix.
#' @details
#' Averaging costs one `n * n` copy per call: repair `d` yourself if calling a
#' solver in a loop.
#' @keywords internal
.Symmetrise <- function(d, tolerance = .SymmetryTolerance()) {
  dev <- SymmetryDeviation_cpp(d)
  if (dev == 0) {
    return(d)
  }
  if (dev > tolerance) {
    stop("`d` must be symmetric: d[i, j] and d[j, i] differ by ",
         format(dev, digits = 3), " relative to their magnitude, beyond the ",
         "tolerance of ", format(tolerance, digits = 3),
         " set by `options(MaxMin.symmetryTolerance=)`")
  }
  warning("`d` is symmetric only to rounding; averaging d[i, j] with d[j, i]. ",
          "Pass `(d + t(d)) / 2` to silence this.", immediate. = TRUE)
  Symmetrised_cpp(d)
}

# Largest scaled discrepancy repaired rather than refused; 0 refuses anything
# inexact. The default is what R's own `isSymmetric()` allows.
.SymmetryTolerance <- function() {
  tol <- getOption("MaxMin.symmetryTolerance", 100 * .Machine$double.eps)
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol < 0) {
    stop("`options(MaxMin.symmetryTolerance=)` must be a single ",
         "non-negative number")
  }
  tol
}

#' Coerce distance input to a square matrix, skipping the round-trip when
#' already a matrix.
#' @param d A `dist` object or a square numeric matrix.
#' @param symmetric Logical: reconcile the two triangles of `d`.
#' @return `.AsDistMatrix()` returns a square numeric matrix.
#' @details
#' A `dist` object is symmetric by construction, so it bypasses the check.
#' @keywords internal
.AsDistMatrix <- function(d, symmetric = TRUE) {
  wasDist <- inherits(d, "dist")
  if (wasDist) {
    d <- as.matrix(d)
  } else if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a `dist` object or a square numeric matrix")
  }
  if (!AllFinite_cpp(d, .NThreads())) {
    stop("distance matrix must not contain NA/NaN/Inf")
  }
  if (symmetric && !wasDist) {
    d <- .Symmetrise(d)
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
#' @return `.MaximinFrom()` returns an integer vector of length `k` of selected row/col indices.
#' @keywords internal
.MaximinFrom <- function(d, k, first) {
  MaximinFrom_cpp(d, as.integer(k), as.integer(first), .NThreads())
}

#' Gonzalez maximin from several starting indices at once
#'
#' Ensemble counterpart of [.MaximinFrom()] / [.MaximinFromPoints()]: solves
#' one greedy pass per seed. Each pass is an independent function of its seed,
#' so under `mc.cores > 1` the passes run concurrently, one per thread; a lone
#' seed, or a problem large enough for a single pass to occupy every thread
#' itself, falls back to the per-pass parallelism instead.
#'
#' @param k Integer: target subsample size (`>= 1`).
#' @param firsts Integer vector of distinct first-selected indices, one per pass.
#' @param d Square pairwise distance matrix, or `NULL` when `points` is given.
#' @param points A `double` `N x dim` coordinate matrix, or `NULL`.
#' @return `.MaximinMulti()` returns a list with one `list(idx, tK)` per element
#'   of `firsts`, in that order.
#' @keywords internal
.MaximinMulti <- function(k, firsts, d = NULL, points = NULL) {
  k <- as.integer(k)
  firsts <- as.integer(firsts)
  res <- if (is.null(points)) {
    MaximinMultiFrom_cpp(d, k, firsts, .NThreads())
  } else {
    MaximinMultiFromPoints_cpp(points, k, firsts, .NThreads())
  }
  lapply(seq_along(firsts), function(j) {
    list(idx = res[["idx"]][, j], tK = res[["t_k"]][[j]])
  })
}

#' Coerce coordinate input for the on-the-fly (matrix-free) samplers
#'
#' The coordinate paths require a complete numeric `N x dim` matrix with
#' `double` storage; the C++ kernels reproduce `stats::dist()`'s exact
#' Euclidean bits, which is only defined for complete data.
#' @param points A numeric matrix (or coercible) of point coordinates.
#' @return `.AsPointsMatrix()` returns a `double` numeric matrix.
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
#' @return `.MaximinFromPoints()` returns an integer vector of selected indices.
#' @keywords internal
.MaximinFromPoints <- function(points, k, first, mask = 0L) {
  MaximinFromPoints_cpp(points, as.integer(k), as.integer(first),
                        as.integer(mask), .NThreads())
}

#' Promote the kernel's free `t_k` score to the user-facing `score` attribute
#'
#' The maximin kernels attach the selection's minimum pairwise distance as a
#' `t_k` attribute (computed during the greedy pass at no extra cost). The
#' ensemble drivers read it via [base::attr()] into `strategy_results`; a bare
#' single pass exposes it directly as the `score` attribute, matching
#' [DropAdd()] and [Grasp()].
#' @param idx Integer vector returned by a maximin kernel.
#' @return `.PromoteScore()` returns `idx` with its `t_k` attribute renamed to `score`.
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
#' @return `.MinPairwiseFromPoints()` returns a numeric scalar; `NA_real_` if `length(idx) < 2`.
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
#' @return `.SubsetScore()` returns a numeric scalar; `NA` if `length(idx) < 2`.
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

#' Greedy farthest-first point selection
#'
#' Greedy farthest-first selection \insertCite{Gonzalez1985,Hochbaum1985}{MaxMin}
#' iteratively selects the point furthest from the current selection to yield a
#' 2-approximation to the _k_-centre and Max Min Diversity problems.
#'
#' @templateVar progress_shows the distance-column path shows a progress bar
#' @template progress
#'
#' @param k Integer: number of points to select.
#' @param d A `dist` object, a square numeric matrix of pairwise distances, or
#' a distance function that takes an index `i` and returns the distance from
#' `i` to each other element (optionally including the self-distance).
#' Ignored when `points` is supplied.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the selection is computed directly from coordinates in
#'   \eqn{O(N \cdot k \cdot dim)} time and \eqn{O(N)} memory, which avoids
#'   creating an \eqn{O(N^2)} distance matrix. Missing entries (`NA`) are not
#'   supported.
#' @param N Integer: the total number of elements. Required (and used) only on
#'   the distance-column oracle path, where it cannot be inferred from the
#'   closure; ignored for the matrix and coordinate paths.
#' @param strategy Integer or character defining how to seed the greedy pass.
#' Pass the name of one or more seeding strategies described in [`PickPoint()`]
#' to run each strategy and return the best solution.
#' @param nSeeds Integer: number of distinct seeds to draw under the (default)
#' `"random_furthest"` strategy. Three starts captures most of the gain from
#' restarting — the improvement curve bends early (knee at n ≈ 3–4 across
#' benchmarks) and additional restarts add little. For higher-quality
#' solutions, prefer [DropAdd()]: tabu search escapes the farthest-first
#' construction family where restarts plateau.
#' @return `FarFirst()` returns an integer vector with class `MaxMinSelection`,
#' listing the selected indices in the order they were selected.
#' Attributes report:
#' - `score`: the selection's minimum pairwise distance (\eqn{T_k}).
#' - `winning_strategy`: character vector listing strategies that attained the
#' optimal score.
#' - `strategy_results`: results for each strategy.
#'
#' @section Parallelism:
#' To parallelize computation when OpenMP is available, set the `"mc.cores"`
#' option:
#'
#' ```r
#' options(mc.cores = 2L)                       # use a fixed number of cores
#' options(mc.cores = parallel::detectCores())  # or all available cores
#' ```
#'
#' The main use case for parallelization is when `nSeeds` is a multiple of
#' `mc.cores`, so each seed point can be evaluated in parallel.
#'
#' @references \insertAllCited{}
#' @seealso [PickPoint()] for the seed indices alone; [DropAdd()] and
#'   [ExactMaxMin()] for higher-effort solvers.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#'
#' # Default: best of three random-furthest starts (set.seed for reproducibility):
#' FarFirst(5L, d)
#'
#' # More random-furthest starts:
#' FarFirst(5L, d, nSeeds = 15L)
#'
#' # Custom two-anchor ensemble:
#' FarFirst(5L, d, strategy = c("diameter", "anti_medoid"))
#'
#' # A single strategy:
#' FarFirst(5L, d, strategy = "diameter")
#'
#' # An explicit start index (integer strategy):
#' FarFirst(5L, d, strategy = 1L)
#'
#' # Matrix-free coordinate path (identical result, O(N) memory):
#' FarFirst(5L, points = pts, strategy = 1L)
#'
#' # Distance-column oracle: supply one column at a time, never the full matrix.
#' data("USArrests")
#' arrestTypes <- USArrests[, c("Murder", "Assault", "Rape")]
#' StateDist <- function(i) {
#'   diffs <- sweep(arrestTypes, 2, unlist(arrestTypes[i, ]), "-")
#'   sqrt(rowSums(diffs ^ 2))
#' }
#' idx <- FarFirst(4L, StateDist, N = nrow(arrestTypes), strategy = 1L)
#' arrestTypes[idx, ]
#' @export
FarFirst <- function(k, d = NULL, points = NULL, N = NULL,
                     strategy = "random_furthest", nSeeds = 3L) {
  progress <- getOption("MaxMin.progress", interactive())
  strategyMissing <- missing(strategy)

  # Every non-empty FarFirst result is a `MaxMinSelection` (a self-describing
  # integer index vector); see `print.MaxMinSelection()`. Stamping it at each
  # exit keeps the matrix, coordinate and oracle paths byte-identical.
  Classify <- function(x) .AsMaxMinSelection(x, "FarFirst")

  # All single-strategy names (the ensemble anchors plus the single-pass-only
  # `"first"`); .kPointEnsembleSeeds is the ensemble-eligible subset.
  validMethods <- c(.kPointEnsembleSeeds, "first")

  if (is.numeric(strategy)) {
    # An integer `strategy` is the explicit 1-based first index (a single bare
    # Gonzalez pass). Guard it here: NA/NaN would otherwise reach the C++
    # kernel as INT_MIN, and a zero-length or multi-element value as an opaque
    # Rcpp "Expecting a single value" error.
    if (length(strategy) != 1L || !is.finite(strategy)) {
      stop("an integer `strategy` (the first index) must be a single finite value")
    }
    first    <- as.integer(strategy)
    strategy <- "first"
  } else if (length(strategy) > 1L) {
    # A multi-element `strategy` requests an ensemble. The downstream
    # match.arg(several.ok = TRUE) silently *drops* names that fail to match,
    # so validate explicitly here against the ensemble anchors (anti_centroid
    # included -- the matrix path drops it later with a warning, not an error).
    bad <- setdiff(strategy, .kPointEnsembleSeeds)
    if (length(bad)) {
      stop("unknown strateg", if (length(bad) > 1L) "ies: " else "y: ",
           paste(bad, collapse = ", "))
    }
    first  <- NULL
  } else {
    first    <- NULL
    strategy <- match.arg(strategy, choices = validMethods)
  }

  # Distance-column oracle path: `d` is a closure returning one matrix column
  # at a time, for metrics with neither a stored matrix nor a coordinate
  # embedding (e.g. on-demand tree-to-tree distances). The selection is
  # identical to the matrix path given the same `first`; only an integer
  # `strategy` (a `first` index) or the deterministic peripheral seed is reachable
  # here, since the richer anchors need O(N^2) work (see Details).
  if (is.function(d)) {
    # A named/character `strategy` is unreachable from an oracle (it would need
    # the whole matrix); warn rather than silently substituting the peripheral
    # seed. `first` is non-NULL only for an integer `strategy`, which *is* honoured.
    if (!strategyMissing && is.null(first)) {
      warning("distance-column oracle path: only an integer `strategy` (a `first` ",
              "index) is honoured; using the deterministic peripheral seed")
    }
    return(Classify(.GonzalezColumn(colFn = d, N = N, k = k, first = first,
                                    progress = progress)))
  }

  if (!is.null(N)) {
    warning("`N` is ignored on the matrix/coordinate path; ",
            "it is only needed when `d` is a distance-column function")
  }
  usePoints <- !is.null(points)
  if (usePoints) {
    points <- .AsPointsMatrix(points)
    nPts <- nrow(points)
  } else {
    # Every Gonzalez path reads whole d(., centre) columns, and the seeding
    # scans cover the full matrix rather than taking a triangle, so d is scored
    # as written and the O(n^2) reconciliation would cost more than the solve.
    d <- .AsDistMatrix(d, symmetric = FALSE)
    nPts <- nrow(d)
  }
  # Pre-check before as.integer(): a finite, length-1, non-negative value.
  # Pre-checking avoids the spurious base-R coercion warning on `k = Inf`.
  if (length(k) != 1L || !is.finite(k) || k < 0L) {
    stop("`k` must be a single non-negative integer")
  }
  k <- min(as.integer(k), nPts)
  if (k == 0L) return(integer(0))

  # The kernels attach a `t_k` attribute (the selection's min pairwise distance,
  # computed for free). A bare single pass exposes it as the `score` attribute
  # (matching DropAdd/Grasp); the ensemble path reads it into strategy_results.
  Greedy <- if (usePoints) {
    function(s) .PromoteScore(.MaximinFromPoints(points, k, first = as.integer(s)))
  } else {
    function(s) .PromoteScore(.MaximinFrom(d, k, first = as.integer(s)))
  }

  # An explicit `first` is a single bare Gonzalez pass, overriding `strategy`.
  if (!is.null(first)) {
    return(Classify(Greedy(first)))
  }

  # The ensemble path runs each named anchor as a full Gonzalez pass and keeps
  # the best by MinDist(). It is taken for a multi-anchor `strategy`, and also for
  # a lone `"random_furthest"` (the default): that token draws `nSeeds` distinct
  # seeds, so it belongs here, not on the single-strategy path. (`PickPoint(
  # strategy = "random_furthest")` still returns exactly one seed for callers who
  # want a single random start.)
  if (length(strategy) > 1L || "random_furthest" %in% strategy) {
    # Matrix path: `"anti_centroid"` is coordinate-only, so drop it (the remaining
    # O(N) seeds cover that role here). Warn only if it was named explicitly,
    # not when filtering the default ensemble.
    anchors <- if (usePoints) strategy else strategy[strategy != "anti_centroid"]
    if (!usePoints && !strategyMissing && length(anchors) < length(strategy)) {
      warning("`anti_centroid` seed requires coordinates; it is dropped on the ",
              "distance-matrix path, where `peripheral` covers the same role")
    }

    nSeeds <- as.integer(nSeeds)
    if (length(nSeeds) != 1L || is.na(nSeeds) || nSeeds < 1L) {
      stop("`nSeeds` must be a single positive integer")
    }

    if (usePoints) {
      return(Classify(.GonzEnsembleFromPoints(points, k, anchors, nSeeds = nSeeds)))
    }
    return(Classify(.GonzEnsemble(d, k, anchors, nSeeds = nSeeds)))
  }

  if (!usePoints && strategy == "anti_centroid") {
    stop("`anti_centroid` seed requires coordinates; supply `points=` or use ",
         "`peripheral` on the distance-matrix path")
  }
  s <- if (usePoints) .PickPoints(points, strategy) else .PickPoint(d, strategy)
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
#' @return `.DistColumn()` returns a numeric vector of length `N`, masked to `-Inf` at position `i`.
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
#' @return `.GonzalezColumn()` returns an integer vector of length `min(k, N)` of selected indices.
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
#' @return `.MaximinFromColumn()` returns an integer vector of selected indices.
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
#' @return `.PeripheralSeedColumn()` returns an integer index of the seed.
#' @keywords internal
.PeripheralSeedColumn <- function(colFn, N) {
  s1 <- which.max(.DistColumn(colFn, 1L, N))
  which.max(.DistColumn(colFn, s1, N))
}
