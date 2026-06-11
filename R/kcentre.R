# kcentre.R
#
# Discrete (vertex) k-centre: choose k centres from the N points so the
# covering radius R(C) = max_p min_{c in C} d(p, c) -- the largest distance from
# any point to its nearest centre -- is minimised. This is the min-MAX covering
# objective of facility location (Gonzalez 1985; Hochbaum & Shmoys 1985),
# distinct from the MMDP/MaxMin objective elsewhere in this package, which
# maximises the min PAIRWISE distance within the selection (see MinDist()).
#
# Three entry points:
#   KCentreRadius()  -- score a centre set by its covering radius (matrix-free
#                       capable, so it scales past the solvers' matrix bound).
#   KCentre()        -- CDSh heuristic (Garcia-Diaz et al. 2017/2019): ~1-3.5%
#                       of optimum at O(n^2 log n), an order of magnitude tighter
#                       than the Gonzalez 2-approximation that FarFirst() gives.
#   ExactKCentre()   -- proven optimum on small instances via a covering integer
#                       program, the covering dual of ExactMaxMin()'s packing IP.

# ----- candidate radii ------------------------------------------------------

# The optimum covering radius is a realised point-to-centre distance, so both
# solvers search the achieved distinct distances. Symmetric `d` is assumed: the
# upper triangle enumerates every off-diagonal value once. The C++ enumerator
# avoids R's n x n `upper.tri` logical mask and `unique`/`order` overhead
# (~28% of a KCentre() call at n = 2000; T-010).
.KCentreCandidates <- function(d) {
  KCentreCandidates_cpp(d)
}

# Deterministic seed vertices for CDSh's i = 0 critical vertex. A peripheral
# Gonzalez pass gives `nstart` well-spread anchors with no RNG; the first is the
# two-sweep peripheral seed.
.KCentrePeripheralSeeds <- function(d, nstart) {
  nstart <- max(1L, as.integer(nstart))
  as.integer(FarFirst(d, min(nstart, nrow(d)), method = "peripheral"))
}

# The k-centre solvers assume a symmetric metric: `KCentreCandidates_cpp` reads
# only the upper triangle and the kernel/covering-IP read d(i,j) as d(j,i).
# `.AsDistMatrix` (unlike the other MaxMin solvers) accepts asymmetric matrices,
# so guard here -- an asymmetric `d` would otherwise give a silently wrong radius
# (KC-002). A `dist` object is symmetric by construction; the O(n^2) check is
# negligible against the O(n^2 log n) solve. (KCentreRadius needs no guard: it
# reads true d(point, centre) columns and is correct for asymmetric input.)
.KCentreRequireSymmetric <- function(d) {
  # In-place C++ check (no transpose copy): negligible beside the O(n^2 log n)
  # solve, where R's isSymmetric() would copy a large dense matrix.
  if (!IsSymmetric_cpp(d, 100 * .Machine$double.eps)) {
    stop("`d` must be symmetric for the k-centre solvers; an asymmetric ",
         "distance would give a silently wrong covering radius")
  }
}

# Above this many candidate radii the exhaustive CDS scan (O(nCand * n^2)) is too
# costly, so KCentre falls back to binary-search CDSh; below it the full scan is
# cheap and avoids the binary search's non-monotone misses (KC-001). The bound
# n(n-1)/2 <= 11325 corresponds to n <= 150.
.kCentreExhaustiveMaxCand <- 11325L

# ----- covering-radius score ------------------------------------------------

#' Covering radius of a centre set (k-centre objective)
#'
#' Returns the covering radius of a set of centres: the largest distance from
#' any of the `N` points to its nearest centre,
#' \eqn{R = \max_p \min_{c \in \mathrm{idx}} d(p, c)}. This is the min-max
#' \emph{k}-centre objective \insertCite{Gonzalez1985}{MaxMin}, the quantity
#' [KCentre()] and [ExactKCentre()] minimise. Lower is better.
#'
#' Unlike [MinDist()] -- the minimum \emph{pairwise} distance \emph{within} a
#' selection (the MMDP objective) -- the covering radius is taken over
#' \emph{all} `N` points and measures how well the centres cover the data.
#'
#' @param d Pairwise distance matrix or `dist` object. Ignored when `points` is
#'   supplied.
#' @param idx Integer vector of centre indices (`>= 1`).
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied the
#'   per-point nearest-centre distances are recomputed from coordinates one
#'   centre column at a time, never materialising the `N x N` matrix (`d` is then
#'   unused), so the score scales to large `N`. For Euclidean data the result is
#'   identical to the matrix path.
#' @return Numeric scalar: the covering radius (`0` when the centres include
#'   every point).
#' @seealso [KCentre()] and [ExactKCentre()] (which minimise this); [MinDist()]
#'   for the complementary MMDP objective.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' centres <- KCentre(d, 4L)
#' KCentreRadius(d, centres)
#' @export
KCentreRadius <- function(d = NULL, idx, points = NULL) {
  idx <- as.integer(idx)
  if (anyNA(idx)) {
    stop("`idx` must not contain NA")
  }
  if (anyDuplicated(idx)) {
    stop("`idx` must not contain duplicate indices")
  }
  if (length(idx) < 1L) {
    stop("`idx` must select at least one centre")
  }
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    if (any(idx < 1L | idx > nrow(points))) {
      stop("`idx` must be centre indices in [1, nrow(points)]")
    }
    nn <- EuclidColFromPoints_cpp(points, idx[1L])
    for (ci in idx[-1L]) {
      nn <- pmin(nn, EuclidColFromPoints_cpp(points, ci))
    }
    return(max(nn))
  }
  d <- .AsDistMatrix(d)
  # Without this guard `d[, 0L]` (or an out-of-range column) yields an empty
  # column and `max()` returns -Inf with only a warning -- the coordinate path
  # already errors. Reject out-of-range indices on both paths (KC-003).
  if (any(idx < 1L | idx > ncol(d))) {
    stop("`idx` must be centre indices in [1, nrow(d)]")
  }
  nn <- d[, idx[1L]]
  for (ci in idx[-1L]) {
    nn <- pmin(nn, d[, ci])
  }
  max(nn)
}

# ----- CDSh heuristic -------------------------------------------------------

#' Near-optimal discrete k-centre by CDSh
#'
#' Chooses `k` centres minimising the covering radius (the largest distance from
#' any point to its nearest centre) with the CDSh heuristic of
#' \insertCite{GarciaDiaz2019;textual}{MaxMin} (see also
#' \insertCite{GarciaDiaz2017}{MaxMin}). CDSh binary-searches the achieved
#' distinct distances; at each trial radius it runs a fixed-`k` farthest-point
#' construction in which every centre is the highest-degree neighbour (within the
#' trial radius) of the currently worst-covered vertex, and accepts the radius
#' when the construction covers all points within it. On the benchmark instances
#' of \insertCite{GarciaDiaz2019;textual}{MaxMin} it reaches roughly 1-3.5% of the
#' optimum at \eqn{O(N^2 \log N)}, far tighter than the Gonzalez 2-approximation
#' that [FarFirst()] gives for this objective (typically tens of per cent above
#' optimum).
#'
#' The achieved covering radius is not monotone in the trial radius, so the binary
#' search can occasionally miss the best candidate. Two safeguards keep the result
#' robust: for a small candidate grid (`n` up to ~150) every radius is scanned
#' exhaustively, and the result is always floored against a deterministic Gonzalez
#' pass, so `KCentre()` is **never worse than the 2-approximation**. For a tighter
#' result at larger `n` raise `nstart`; for the proven optimum on a small instance
#' use [ExactKCentre()].
#'
#' The construction is otherwise fully deterministic: where the reference seeds its
#' first critical vertex at random, `KCentre()` uses deterministic peripheral
#' anchors (`nstart` of them, the best kept). Like [ExactKCentre()] this is a
#' distance-matrix method, \eqn{O(N^2)} in memory; for the covering radius of an
#' existing selection at larger `N`, [KCentreRadius()] has a matrix-free path.
#'
#' @param d A `dist` object or a square symmetric numeric distance matrix.
#' @param k Integer number of centres, `1 <= k <= nrow(d)`.
#' @param nstart Integer; how many deterministic peripheral seeds to try, keeping
#'   the lowest-radius result. Default `1`. Ignored if `seeds` is supplied.
#' @param seeds Optional integer vector of explicit 1-based seed vertices for the
#'   construction's first critical vertex (overrides `nstart`).
#' @return Integer vector of length `<= k` (ascending): the chosen centres. The
#'   achieved covering radius is attached as attribute `radius`. The vector has
#'   class `"KCentreSelection"` and prints as a one-line summary; it is otherwise
#'   an ordinary integer vector and indexes a matrix or coordinate set directly.
#' @seealso [ExactKCentre()] for the proven optimum; [KCentreRadius()] to score a
#'   selection; [FarFirst()] for the Gonzalez 2-approximation baseline.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(120), ncol = 2)
#' d <- dist(pts)
#' centres <- KCentre(d, 5L)
#' KCentreRadius(d, centres)
#' # CDSh covers at least as tightly as the Gonzalez 2-approximation:
#' KCentreRadius(d, FarFirst(d, 5L, method = "peripheral"))
#' @export
KCentre <- function(d, k, nstart = 1L, seeds = NULL) {
  needSymCheck <- !inherits(d, "dist")          # a dist object is symmetric
  d <- .AsDistMatrix(d)
  if (needSymCheck) .KCentreRequireSymmetric(d)
  n <- nrow(d)
  if (length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("`k` must be a single positive integer")
  }
  k <- as.integer(k)
  if (k > n) {
    stop("`k` must satisfy 1 <= k <= nrow(d)")
  }
  # k == n: every point is its own centre; the covering radius is 0.
  if (k == n) {
    return(structure(seq_len(n), radius = 0, producer = "CDSh",
                     class = "KCentreSelection"))
  }
  cand <- .KCentreCandidates(d)
  if (is.null(seeds)) {
    seeds <- .KCentrePeripheralSeeds(d, nstart)
  } else {
    seeds <- as.integer(seeds)
    if (anyNA(seeds) || any(seeds < 1L | seeds > n)) {
      stop("`seeds` must be indices in [1, nrow(d)]")
    }
  }
  # The per-radius achieved covering radius is non-monotone in r, so the O(log n)
  # binary search can skip the best candidate. When the candidate grid is small
  # enough, scan it exhaustively instead (KC-001).
  exhaustive <- length(cand) <= .kCentreExhaustiveMaxCand
  best <- NULL
  bestR <- Inf
  for (s in seeds) {
    res <- KCentreCDSh_cpp(d, k, s, cand, exhaustive)
    r <- attr(res, "radius")
    if (r < bestR) {
      bestR <- r
      best <- as.integer(res)
    }
  }
  # Gonzalez floor (KC-001): the binary-search CDSh can occasionally return a
  # covering radius ABOVE the Gonzalez 2-approximation it is documented to beat
  # (the feasibility predicate is non-monotone). Take the better of CDSh and a
  # cheap deterministic Gonzalez pass, guaranteeing the result is never worse
  # than the 2-approximation at O(n*k) extra cost.
  gonz <- as.integer(FarFirst(d, k, method = "peripheral"))
  gonzR <- KCentreRadius(d, gonz)
  if (gonzR < bestR) {
    bestR <- gonzR
    best <- gonz
  }
  # The construction can occasionally re-pick an already-chosen vertex as a
  # centre; a duplicate never changes the covering radius, so collapse to the
  # distinct centre set (which may then be smaller than k).
  structure(sort(unique(best)), radius = bestR, producer = "CDSh",
            class = "KCentreSelection")
}

# ----- exact covering IP ----------------------------------------------------

# Solve one minimum-cover feasibility probe at radius r and classify it against
# the centre budget k. The IP minimises the number of open centres subject to
# every point being within r of an open centre; r is feasible for k-centre iff
# that minimum is <= k. The witness is validated independently of the solver
# status (the chosen centres are checked to cover every point within r), exactly
# as .MaxISVerdict() validates its independent set. Verdicts mirror that helper:
# "feasible" (a <= k cover found), "infeasible" (IP proven optimal, min cover
# > k), or "inconclusive" (budget expired).
.MinCoverVerdict <- function(d, n, r, k, timeLimit) {
  if (!is.finite(timeLimit) || timeLimit <= 0) {
    return(list(verdict = "inconclusive", witness = integer(0)))
  }
  # Coverage incidence at radius r: centre i covers point j iff d(i, j) <= r
  # (includes i = j, distance 0). `which(arr.ind)` columns are (row = centre i,
  # col = point j); the covering constraint for point j sums y over its centres.
  cover <- which(d <= r, arr.ind = TRUE)
  A <- Matrix::sparseMatrix(i = cover[, 2L], j = cover[, 1L], x = 1,
                            dims = c(n, n))
  res <- highs::highs_solve(
    L       = rep.int(1, n),
    lower   = rep.int(0, n),
    upper   = rep.int(1, n),
    A       = A,
    lhs     = rep.int(1, n),                  # each point covered at least once
    rhs     = rep.int(Inf, n),
    types   = rep.int("I", n),                # integer var on [0, 1] = binary
    maximum = FALSE,
    control = list(
      threads    = 1L,                        # determinism
      time_limit = timeLimit
    )
  )
  sel <- which(res$primal_solution > 0.5)
  # Independent validation: do the chosen centres cover every point within r?
  validCover <- length(sel) >= 1L &&
    all(apply(d[sel, , drop = FALSE] <= r, 2L, any))

  if (validCover && length(sel) <= k) {
    return(list(verdict = "feasible", witness = sel))
  }
  # Infeasibility is provable only when the IP reached optimality and the
  # certified minimum cover still exceeds k.
  optimal <- identical(res$status_message, "Optimal")
  if (optimal && validCover && length(sel) > k) {
    return(list(verdict = "infeasible", witness = integer(0)))
  }
  list(verdict = "inconclusive", witness = integer(0))
}

#' Exact discrete k-centre optimum on small instances
#'
#' Solves the discrete (vertex) \emph{k}-centre problem to proven optimality: the
#' optimum covering radius is the smallest threshold `r`, over the achieved
#' distinct distances, for which `k` centres can cover every point within `r`.
#' Each probe solves a minimum-cardinality \emph{set-cover} integer program with
#' the `highs` MILP backend, the covering constraints held as a sparse matrix --
#' the covering dual of [ExactMaxMin()]'s node-packing program. The search is
#' warm-started from the [KCentre()] (CDSh) radius, a proven feasible upper bound
#' that caps the binary search, then bisects downward to the smallest feasible
#' radius.
#'
#' The covering optimum may be attained by fewer than `k` centres (extra centres
#' never help once coverage is achieved); `indices` then has length `< k` and the
#' reported `radius` is still the proven \emph{k}-centre optimum. The problem is
#' NP-hard, so this is an external ground-truth reference for small instances,
#' not a scalable method.
#'
#' @param d A `dist` object or a square symmetric numeric distance matrix.
#' @param k Integer centre budget, `1 <= k <= nrow(d)`.
#' @param solver Solver to use. Currently only `"highs"` is implemented; `NULL`
#'   selects it.
#' @param maxSeconds Wall-clock budget in seconds for the whole search (shared
#'   across the internal IP solves). If it expires before the optimum is proven,
#'   the smallest radius proven feasible so far is returned with `proven = FALSE`.
#' @param warmStart Currently unused; reserved for a caller-supplied feasible
#'   centre set. The internal CDSh warm start runs regardless.
#' @param progress Logical; show a progress indicator. Default: `TRUE` in
#'   interactive sessions (`getOption("MaxMin.progress", interactive())`).
#' @return A list of class `"KCentreExact"` with fields
#'   \describe{
#'     \item{indices}{Integer vector (ascending), length `<= k`: the centres.}
#'     \item{radius}{The covering radius they achieve; the proven optimum when
#'       `proven` is `TRUE`, otherwise a valid upper bound.}
#'     \item{proven}{Logical: `TRUE` if optimality was certified within budget.}
#'     \item{time_s}{Wall-clock seconds elapsed.}
#'     \item{solver}{Name of the MILP backend used.}
#'     \item{n, k}{Instance size and centre budget.}
#'     \item{n_centres}{`length(indices)`.}
#'   }
#'   It prints as a one-line summary and is otherwise an ordinary list.
#' @seealso [KCentre()] for the fast near-optimal heuristic; [ExactMaxMin()] for
#'   the dual MMDP optimum.
#' @references \insertAllCited{}
#' @examples
#' \donttest{
#' if (requireNamespace("highs", quietly = TRUE) &&
#'     requireNamespace("Matrix", quietly = TRUE)) {
#'   set.seed(1)
#'   pts <- matrix(rnorm(40), ncol = 2)
#'   d <- dist(pts)
#'   ExactKCentre(d, 3L)
#' }
#' }
#' @export
ExactKCentre <- function(d, k, solver = NULL, maxSeconds = 60,
                         warmStart = NULL,
                         progress = getOption("MaxMin.progress", interactive())) {
  t0 <- proc.time()[[3L]]
  if (is.null(solver)) solver <- "highs"
  if (!identical(solver, "highs")) {
    stop("Unsupported `solver`: ", solver, ". Only \"highs\" is implemented.")
  }
  if (!requireNamespace("highs", quietly = TRUE)) {
    stop("The `highs` package is required for ExactKCentre(). ",
         "Install it with install.packages(\"highs\").")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The `Matrix` package is required for ExactKCentre(). ",
         "Install it with install.packages(\"Matrix\").")
  }

  needSymCheck <- !inherits(d, "dist")          # a dist object is symmetric
  d <- .AsDistMatrix(d)
  if (needSymCheck) .KCentreRequireSymmetric(d)
  n <- nrow(d)
  if (length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("`k` must be a single positive integer")
  }
  k <- as.integer(k)
  if (k > n) {
    stop("`k` must satisfy 1 <= k <= nrow(d)")
  }
  Elapsed <- function() proc.time()[[3L]] - t0

  # `radius` is always the true covering radius of the returned centres, computed
  # from the witness rather than the search threshold (KC-004): for a proven
  # optimum it equals the smallest feasible candidate, and for an unproven /
  # fallback witness it is the tight achieved radius, so `res$radius` and
  # `KCentreRadius(d, res$indices)` never disagree.
  Pack <- function(indices, proven) {
    indices <- sort(as.integer(indices))
    structure(
      list(indices = indices, radius = KCentreRadius(d, indices), proven = proven,
           time_s = Elapsed(), solver = solver, n = n, k = k,
           n_centres = length(indices)),
      class = "KCentreExact"
    )
  }

  # k == n: every point is its own centre; covering radius 0, trivially proven.
  if (k == n) {
    return(Pack(seq_len(n), TRUE))
  }

  cand <- .KCentreCandidates(d)
  feas <- function(idx, remaining) .MinCoverVerdict(d, n, cand[idx], k, remaining)

  if (progress) {
    .pb <- cli::cli_progress_bar("ExactKCentre", total = NA, .auto_close = FALSE)
  }
  tick <- function() if (progress) cli::cli_progress_update(id = .pb)

  # Warm start: CDSh gives k centres with a feasible covering radius, so the
  # optimum index is at or below that radius. The CDSh witness is the incumbent.
  ws <- tryCatch(KCentre(d, k), error = function(e) NULL)
  if (is.null(ws)) {
    # Fallback: any single point covers all within the diameter (the largest
    # candidate), a loose but valid feasible upper bound.
    hi <- length(cand); bestIdx <- hi; bestWitness <- 1L
  } else {
    hi <- findInterval(attr(ws, "radius"), cand)   # cand[hi] == achieved radius
    bestIdx <- hi; bestWitness <- as.integer(ws)
  }

  # Smallest-feasible bisection over (0, hi]: feasibility is monotone (a larger
  # radius only enlarges every point's covering set, so the minimum cover can
  # only shrink), so the feasible candidates form a suffix and cand[hi] is known
  # feasible from the warm start.
  inconclusive <- FALSE
  lo <- 1L
  while (lo < hi) {
    rem <- maxSeconds - Elapsed()
    if (rem <= 0) { inconclusive <- TRUE; break }
    mid <- (lo + hi) %/% 2L
    v <- feas(mid, rem); tick()
    if (identical(v$verdict, "feasible")) {
      hi <- mid; bestIdx <- mid; bestWitness <- v$witness
    } else if (identical(v$verdict, "infeasible")) {
      lo <- mid + 1L
    } else {
      inconclusive <- TRUE; break
    }
  }

  if (progress) cli::cli_progress_done(id = .pb)

  Pack(bestWitness, !inconclusive)
}

# ----- S3 display -----------------------------------------------------------

#' Format and print k-centre solver results
#'
#' One-line summaries of the objects returned by [KCentre()]
#' (`"KCentreSelection"`) and [ExactKCentre()] (`"KCentreExact"`): the centre
#' count, the chosen indices, the method, and the achieved covering radius (with
#' proof status for the exact solver). Both objects are otherwise unchanged.
#' @param x A `"KCentreSelection"` or `"KCentreExact"` object.
#' @param ... Ignored; present for S3 compatibility.
#' @return `x`, invisibly (`print`); a length-1 character string (`format`).
#' @name print.KCentre
#' @family reporting functions
#' @examples
#' set.seed(1)
#' KCentre(dist(matrix(rnorm(60), ncol = 2)), 4L)
#' @export
format.KCentreSelection <- function(x, ...) {
  idx <- as.integer(x)
  nc <- length(idx)
  sprintf("%d centre%s (%s) by CDSh, covering radius <= %s",
          nc, if (nc == 1L) "" else "s", .FormatIndexList(idx),
          format(signif(attr(x, "radius"), 4L)))
}

#' @rdname print.KCentre
#' @export
print.KCentreSelection <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  invisible(x)
}

#' @rdname print.KCentre
#' @export
format.KCentreExact <- function(x, ...) {
  nc <- length(x$indices)
  status <- if (isTRUE(x$proven)) {
    sprintf("exact MILP (%s), proven optimal", x$solver)
  } else {
    sprintf("exact MILP (%s), unproven incumbent", x$solver)
  }
  rel <- if (isTRUE(x$proven)) "=" else "<="
  sprintf("%d centre%s (%s) by %s, covering radius %s %s",
          nc, if (nc == 1L) "" else "s", .FormatIndexList(x$indices),
          status, rel, format(signif(x$radius, 4L)))
}

#' @rdname print.KCentre
#' @export
print.KCentreExact <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  invisible(x)
}

# ----- US-spelling aliases --------------------------------------------------
# UK spelling (KCentre / ExactKCentre / KCentreRadius) is canonical; these
# expose the US spelling of each exported k-centre function as an identical
# alias, so callers using either spelling reach the same function.

#' @rdname KCentreRadius
#' @export
KCenterRadius <- KCentreRadius

#' @rdname KCentre
#' @export
KCenter <- KCentre

#' @rdname ExactKCentre
#' @export
ExactKCenter <- ExactKCentre
