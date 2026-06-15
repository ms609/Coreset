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
  as.integer(FarFirst(min(nstart, nrow(d)), d, strategy = "peripheral"))
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

#' Covering radius of a set of centres
#'
#' `KCentreRadius()` computes the covering radius of a set of centres:
#' the largest distance from any of the `N` points to its nearest centre,
#' \eqn{R = \max_p \min_{c \in \mathrm{idx}} d(p, c)}. This is the min-max
#' \emph{k}-centre objective \insertCite{Gonzalez1985}{MaxMin} minimized by
#' [KCentre()] and [ExactKCentre()].
#'
#' @param d Pairwise distance matrix or `dist` object. Ignored when `points` is
#'   supplied.
#' @param idx Integer vector of centre indices (`>= 1`).
#' @param points `N x dim` numeric coordinate matrix. When supplied, the
#'  per-point nearest-centre distances are computed from coordinates one
#'  column at a time. This supports sets with large `N`, in which the full
#'  `N x N` matrix would overflow available memory.
#' @return `KCentreRadius()` returns a numeric denoting the covering radius.
#' @seealso [KCentre()] and [ExactKCentre()] (which minimise this); [MinDist()]
#'   for the complementary MMDP objective.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' centres <- KCentre(4L, d)
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

#' Near-optimal discrete k-centre solver
#'
#' `KCentre()` selects $k$ elements (centres) so as to minimize the largest
#' distance from any point to its nearest centre (the covering radius),
#' using the CDSh heuristic
#' \insertCite{GarciaDiaz2017,GarciaDiaz2019}{MaxMin}.
#'
#' On the benchmark instances of \insertCite{GarciaDiaz2019;textual}{MaxMin},
#' CDSh reaches roughly 1-3.5% of the optimum at \eqn{O(N^2 \log N)},
#' far tighter than [FarFirst()] (typically tens of per cent above optimum).
#'
#' The achieved covering radius is not monotone in the trial radius, so the CDSh
#' binary search can occasionally miss the best candidate.
#' As a safeguard, `KCentre()` runs by default an exhaustive search of
#' a small candidate grid (for `n` up to ~150); and an additional `FarFirst()`
#' pass (controlled via the `effort` argument). These safeguards ensure
#' that `KCentre()` always returns at least a 2-approximation.
#'
#' @param k Integer specifying maximum number of centres to identify,
#' from 1 to `nrow(d)`.
#' @param d `dist` object or a square symmetric numeric distance matrix.
#' @param nstart Integer specifying how many deterministic peripheral seeds to
#'  try.
#' @param effort Integer: if `> 0`, run a parallel `FarFirst()` search
#'  with `effort` random seeds, returning the best of all results.
#'
#' @return `KCentre()` returns an integer vector of length `<= k` specifying
#' the chosen centres in ascending order.
#' The achieved covering radius is attached as attribute `radius`.
#' The vector has class `"KCentreSelection"` and prints as a one-line summary.
#' @seealso [ExactKCentre()] for the proven optimum;
#' [KCentreRadius()] for a selection's score;
#' [FarFirst()] for the \insertCite{Gonzalez1985;textual}{MaxMin}
#' 2-approximation baseline.
#' @references \insertAllCited{}
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(120), ncol = 2)
#' d <- dist(pts)
#' centres <- KCentre(5L, d)
#' KCentreRadius(d, centres)
#' # Results will beat a Gonzalez 2-approximation:
#' KCentreRadius(d, FarFirst(5L, d, nSeeds = 1))
#' @export
KCentre <- function(k, d, nstart = 1L, effort = 1L) {
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
  effort <- as.integer(effort)
  if (length(effort) != 1L || is.na(effort) || effort < 0L) {
    stop("`effort` must be a single non-negative integer")
  }
  # k == n: every point is its own centre; the covering radius is 0.
  if (k == n) {
    return(structure(seq_len(n), radius = 0, producer = "CDSh",
                     class = "KCentreSelection"))
  }
  cand <- .KCentreCandidates(d)
  seeds <- .KCentrePeripheralSeeds(d, nstart)
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
  # Gonzalez floor (KC-001), controlled by `effort`. The binary-search CDSh can
  # occasionally return a covering radius ABOVE the Gonzalez 2-approximation it is
  # documented to beat (the feasibility predicate is non-monotone), so by default
  # the result is floored against a Gonzalez pass -- never worse than the
  # 2-approximation. `effort = 0` skips the floor (raw CDSh, no guarantee);
  # `effort = 1` runs one deterministic peripheral pass; `effort > 1` runs a
  # distinct-seed restart with `nSeeds = effort` (a tighter floor, drawing on the
  # session RNG).
  if (effort >= 1L) {
    gonz <- if (effort == 1L) {
      as.integer(FarFirst(k, d, strategy = "peripheral"))
    } else {
      as.integer(FarFirst(k, d, nSeeds = effort))
    }
    gonzR <- KCentreRadius(d, gonz)
    if (gonzR < bestR) {
      bestR <- gonzR
      best <- gonz
    }
  }
  # Centres are distinct by construction (an already-chosen vertex is never
  # re-picked); unique() collapses only the degenerate fully-covered case
  # (radius 0), where fewer than k distinct centres already cover everything.
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
  if (!is.finite(timeLimit) || timeLimit <= 0) { # nocov start
    return(list(verdict = "inconclusive", witness = integer(0)))
  } # nocov end
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
  list(verdict = "inconclusive", witness = integer(0))  # nocov
}

#' Exact discrete k-centre optimum on small instances
#'
#' `ExactKCentre()` finds an optimal solution to the discrete \emph{k}-centre
#' problem.
#'
#' The optimum covering radius is the smallest threshold `r`, over the achieved
#' distinct distances, for which `k` centres can cover every point within `r`.
#'
#' Each probe solves a minimum-cardinality \emph{set-cover} integer program with
#' the `highs` MILP backend, the covering constraints held as a sparse matrix --
#' the covering dual of [ExactMaxMin()]'s node-packing program. The search is
#' warm-started from the [KCentre()] (CDSh) radius, a proven feasible upper bound
#' that caps the binary search, then bisects downward to the smallest feasible
#' radius.
#'
#' @param k Integer specifying maximum number of centres, `1 <= k <= nrow(d)`.
#' @param d A `dist` object or a square symmetric numeric distance matrix.
#' @param maxSeconds Wall-clock budget in seconds for the whole search.
#' If it expires before the optimum is proven, the smallest radius proven
#' feasible so far is returned, with the attribute `proven = FALSE`.
#' @templateVar progress_shows a progress indicator is shown
#' @template progress
#' @return `ExactKCentre()` returns an integer vector of length `<= k`
#'   (ascending): the chosen centres. It has class
#'   `c("KCentreExact", "KCentreSelection")` and the following attributes:
#'   \describe{
#'     \item{radius}{The covering radius they achieve; the proven optimum when
#'       `proven` is `TRUE`, otherwise a valid upper bound.}
#'     \item{proven}{Logical: `TRUE` if optimality was certified within budget.}
#'     \item{time_s}{Wall-clock seconds elapsed.}
#'     \item{n, k}{Instance size and centre budget.}
#'     \item{n_centres}{`length(result)`.}
#'   }
#'   It prints as a one-line summary and indexes a matrix or data frame directly.
#'   The `"KCentreSelection"` superclass means [KCentreRadius()] and any generic
#'   written for that class work here too.
#'
#' The covering optimum may be attained by fewer than `k` centres (extra centres
#' never help once coverage is achieved); the result then has length `< k` and the
#' reported `radius` is still the proven \emph{k}-centre optimum. The problem is
#' NP-hard, so this is an external ground-truth reference for small instances,
#' not a scalable method.
#'
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
#'   ExactKCentre(3L, d)
#' }
#' }
#' @export
ExactKCentre <- function(k, d, maxSeconds = 60) {
  progress <- getOption("MaxMin.progress", interactive())
  t0 <- proc.time()[[3L]]
  if (!requireNamespace("highs", quietly = TRUE)) { # nocov start
    stop("The `highs` package is required for ExactKCentre(). ",
         "Install it with install.packages(\"highs\").")
  } # nocov end
  if (!requireNamespace("Matrix", quietly = TRUE)) { # nocov start
    stop("The `Matrix` package is required for ExactKCentre(). ",
         "Install it with install.packages(\"Matrix\").")
  } # nocov end

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
  # fallback witness it is the tight achieved radius, so `attr(res, "radius")`
  # and `KCentreRadius(d, as.integer(res))` never disagree.
  Pack <- function(indices, proven) {
    indices <- sort(as.integer(indices))
    structure(
      indices,
      radius    = KCentreRadius(d, indices),
      proven    = proven,
      time_s    = Elapsed(),
      solver    = "highs",
      n         = n,
      k         = as.integer(k),
      n_centres = length(indices),
      class     = c("KCentreExact", "KCentreSelection")
    )
  }

  # k == n: every point is its own centre; covering radius 0, trivially proven.
  if (k == n) {
    return(Pack(seq_len(n), TRUE))
  }

  cand <- .KCentreCandidates(d)
  feas <- function(idx, remaining) .MinCoverVerdict(d, n, cand[idx], k, remaining)

  if (progress) {
    .pb <- cli::cli_progress_bar("ExactKCentre", total = NA, .auto_close = FALSE) # nocov
  }
  tick <- function() if (progress) cli::cli_progress_update(id = .pb)

  # Warm start: CDSh gives k centres with a feasible covering radius, so the
  # optimum index is at or below that radius. The CDSh witness is the incumbent.
  ws <- tryCatch(KCentre(k, d), error = function(e) NULL)
  if (is.null(ws)) { # nocov start
    # Fallback: any single point covers all within the diameter (the largest
    # candidate), a loose but valid feasible upper bound.
    hi <- length(cand); bestIdx <- hi; bestWitness <- 1L
  } else { # nocov end
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
    if (rem <= 0) { inconclusive <- TRUE; break } # nocov
    mid <- (lo + hi) %/% 2L
    v <- feas(mid, rem); tick()
    if (identical(v$verdict, "feasible")) {
      hi <- mid; bestIdx <- mid; bestWitness <- v$witness
    } else if (identical(v$verdict, "infeasible")) {
      lo <- mid + 1L
    } else {  # nocov start
      inconclusive <- TRUE; break
    }  # nocov end
  }

  if (progress) cli::cli_progress_done(id = .pb) # nocov

  Pack(bestWitness, !inconclusive)
}

# ----- US spelling aliases --------------------------------------------------

#' @rdname KCentreRadius
#' @export
KCenterRadius <- KCentreRadius

#' @rdname KCentre
#' @export
KCenter <- KCentre

#' @rdname ExactKCentre
#' @export
ExactKCenter <- ExactKCentre
