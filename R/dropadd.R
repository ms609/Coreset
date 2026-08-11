# dropadd.R
#
# DropAdd Tabu Search (DropAdd-TS) for the Max-Min Diversity Problem.
#
# Implementation of Porumbel, Hao & Glover (2011) "A simple and effective
# algorithm for the MaxMin diversity problem", Annals of Operations Research
# 186:275-293. Algorithms 1-4 of that paper.
#
# The MMDPo (eq. 2 of the paper) objective is
#     max  [ min_{x,y in X} d(x,y) + eps * sum_{x,y in X} d(x,y) ]
# with a small eps. We use eps = 1e-9. The first term is the canonical MMDP
# (MaxMin) objective; the sum acts purely as a tie-break.
#
# Used in the manuscript as the SOTA heuristic competitor for MMDP per Marti
# et al. 2022 EJOR.

# ----- helpers --------------------------------------------------------------

# Constructive phase (Algorithm 1).
# Returns list(S = integer(k), iter_add = integer(k)) where iter_add[i] is
# the iteration at which S[i] was added (= i).
.DropAddConstruct <- function(dmat, k) {
  n <- nrow(dmat)
  S <- integer(k)
  # Seed point: argmax over Z of sum_y d(x, y) (Porumbel's eq. before Alg. 1).
  rowSumsD <- rowSums(dmat)
  S[1L] <- which.max(rowSumsD)  # ties: which.max → smallest index (deterministic)
  inS <- logical(n)
  inS[S[1L]] <- TRUE

  # Streamlined records.
  # minDist(x):       min over y in current_S of d(x, y), excluding x itself
  # sumDist(x):       sum over y in current_S of d(x, y) (d(x,x)=0 so harmless)
  # minDistCount(x): | { y in current_S, y != x : d(x, y) == minDist(x) } |
  minDist       <- dmat[, S[1L]]            # length n; minDist(S[1])=0 (self)
  minDist[S[1L]] <- Inf                     # never want to compare a selected
                                             # point to itself in min calcs
  sumDist       <- dmat[, S[1L]]
  minDistCount <- rep(1L, n)
  minDistCount[S[1L]] <- 0L                # no "other selected" yet

  if (k >= 2L) {
    for (h in 2L:k) {
      # Among AddX = Z \ S, pick x maximising (minDist, sumDist) lex,
      # ties → smallest index.
      cand <- which(!inS)
      md <- minDist[cand]
      bestMin <- max(md)
      tied <- cand[md == bestMin]
      if (length(tied) > 1L) {
        sd <- sumDist[tied]
        bestSum <- max(sd)
        tied <- tied[sd == bestSum]
      }
      xNew <- tied[1L]   # smallest index among lex-ties
      S[h] <- xNew
      inS[xNew] <- TRUE

      # Update streamlined records for ADD (Algorithm 3).
      # For every x != xNew:
      dCol <- dmat[, xNew]
      others <- seq_len(n)[-xNew]
      dToNew <- dCol[others]
      # sumDist
      sumDist[others] <- sumDist[others] + dToNew
      # Cases
      mdO <- minDist[others]
      caseB <- dToNew < mdO
      caseA <- dToNew == mdO & !caseB
      # Apply
      if (any(caseB)) {
        idx <- others[caseB]
        minDist[idx] <- dToNew[caseB]
        minDistCount[idx] <- 1L
      }
      if (any(caseA)) {
        idx <- others[caseA]
        minDistCount[idx] <- minDistCount[idx] + 1L
      }
      # xNew itself: its minDist becomes min over (S - {xNew}) of d.
      # Since we'd already maintained minDist[xNew] = Inf, we need to set
      # it now to the min distance from xNew to the previously-selected
      # points (which is what its minDist would have been without the
      # Inf override).
      prevS <- S[seq_len(h - 1L)]
      dXnew <- dmat[xNew, prevS]
      mn <- min(dXnew)
      minDist[xNew] <- mn
      minDistCount[xNew] <- sum(dXnew == mn)
    }
  }

  list(S = S, inS = inS, minDist = minDist, sumDist = sumDist,
       minDistCount = minDistCount)
}

# ----- distance-column oracle path ------------------------------------------
#
# The oracle path is the pure-R twin of src/dropadd.cpp: identical streamlined
# records (minDist / sumDist / minDistCount), identical FIFO circular-buffer
# drop, identical lexicographic (minDist, sumDist) add with the just-dropped
# x# excluded for one iteration, identical MMDPo objective. The ONLY change is
# the distance source -- each needed column d(., x) comes from a caller-supplied
# `colFn(x)` instead of a matrix read `dmat[, x]`, so the N x N matrix is never
# materialised. This is the same "same algorithm, different distance source"
# framing as src/dropadd_mf.cpp, and the counterpart for DropAdd of
# .MaximinFromColumn() for FarFirst.
#
# Oracle-call budget: k for the construction (the seed's own column, then one
# per added element), plus 2 for the default peripheral seed, plus 2 per
# main-loop iteration, plus one per entry of `needRecompute`. Unlike the matrix
# path the per-iteration cost is *distance evaluations*, so `plateau` and
# `maxSeconds` control the oracle bill directly; see the DropAdd() roxygen for
# the caller-side closure advice.
#
# Symmetry: as on FarFirst()'s oracle path, the records assume d(i, j) = d(j, i).
# The matrix path makes the same assumption (.AsDistMatrix() deliberately skips
# the O(N^2) symmetry check).

#' Normalise a distance-column oracle result for the DropAdd records
#'
#' [.DistColumn()] masks the self-distance to `-Inf`, which is what
#' `which.max()` / `pmin.int()` want on [FarFirst()]'s oracle path. DropAdd's
#' streamlined records want the opposite: `src/dropadd.cpp` reads a matrix whose
#' diagonal is `0` and folds the *whole* column into `sumDist` (relying on
#' `d(x, x) = 0` being harmless), masking `minDist` at the selected index
#' separately. A `-Inf` self entry would poison `sumDist` irrecoverably and
#' corrupt the MMDPo tie-break while leaving the returned indices plausible. So
#' take [.DistColumn()]'s length normalisation and `NA` guard, then restore the
#' zero diagonal.
#' @inheritParams .DistColumn
#' @return `.DropAddColumn()` returns a numeric vector of length `N` whose
#'   position `i` is `0`, matching a distance matrix's diagonal.
#' @keywords internal
.DropAddColumn <- function(colFn, i, N) {
  col <- .DistColumn(colFn, i, N)
  col[i] <- 0                        # restore the matrix diagonal
  # Return:
  col
}

#' Choose the next element to add to a DropAdd selection
#'
#' The `argmax` over the unselected elements of `(minDist, sumDist)` taken
#' lexicographically, ties broken to the smallest index. This is exactly the
#' rule of the add scans in `src/dropadd.cpp`: a single ascending pass that
#' replaces the incumbent only on a strict improvement, so it lands on the first
#' index attaining the lexicographic maximum.
#' @param st Record environment; see [.DropAddConstructColumn()].
#' @param exclude Integer index barred from selection this iteration (`0L` for
#'   none). Used for the just-dropped `x#`, which Porumbel et al. (2011, p. 281)
#'   exclude from `Add X(k)` for one iteration -- the tabu rule.
#' @return `.DropAddPick()` returns a single integer index.
#' @keywords internal
.DropAddPick <- function(st, exclude = 0L) {
  eligible <- !st$inS
  if (exclude > 0L) {
    eligible[exclude] <- FALSE
  }
  cand <- which(eligible)
  md <- st$minDist[cand]
  tied <- cand[md == max(md)]
  if (length(tied) > 1L) {
    sd <- st$sumDist[tied]
    tied <- tied[sd == max(sd)]
  }
  # Return:
  tied[[1L]]
}

#' Fold a newly added element into the streamlined DropAdd records
#'
#' The ADD pass of Algorithm 3, shared by the construction and the tabu loop and
#' mirroring the add-pass blocks of `src/dropadd.cpp`. Mutates `st` in place
#' (the package forbids `<<-`; `st` is a `new.env(parent = emptyenv())` record
#' bundle). `xNew`'s *own* record is set by the caller, which knows which peers
#' are current.
#' @inheritParams .DropAddPick
#' @param col Self-zeroed distance column of `xNew`; see [.DropAddColumn()].
#' @param xNew Integer index of the element just added.
#' @return `.DropAddApplyAdd()` returns `NULL` invisibly, for its effect on `st`.
#' @keywords internal
.DropAddApplyAdd <- function(st, col, xNew) {
  st$sumDist <- st$sumDist + col       # d(xNew, xNew) = 0, harmless
  # Both masks are read off the pre-update minDist, and are disjoint, so the
  # vectorised form matches the C++'s in-place ascending scan exactly.
  lower <- col < st$minDist
  equal <- col == st$minDist
  lower[xNew] <- FALSE
  equal[xNew] <- FALSE
  st$minDist[lower] <- col[lower]
  st$minDistCount[lower] <- 1L
  st$minDistCount[equal] <- st$minDistCount[equal] + 1L
  # Return:
  invisible(NULL)
}

#' Constructive phase of DropAdd from a distance-column oracle
#'
#' Algorithm 1 of \insertCite{Porumbel2011;textual}{MaxMin} against a column
#' oracle: the counterpart of `.DropAddConstruct()` (and of the construction
#' block of `src/dropadd.cpp`) that never materialises the distance matrix.
#' Costs one oracle call per selected element.
#'
#' @section Seed deviation:
#' The matrix kernel seeds at the max-row-sum point,
#' \eqn{\mathrm{argmax}_x \sum_y d(x, y)}. From a column oracle that needs all
#' `N` columns -- the full \eqn{O(N^2)} pairwise set, the very cost this path
#' exists to avoid. We therefore substitute the \eqn{O(N)} two-sweep peripheral
#' point of [.PeripheralSeedColumn()] (the element furthest from element 1, then
#' the element furthest from that; two oracle calls, no RNG). This is the
#' column-oracle twin of the centroid-peripheral seed substituted on the
#' coordinate path (`src/dropadd_mf.cpp`), and rests on the same argument: a
#' peripheral element has a large total distance to the rest, so it tracks the
#' max-row-sum point closely. The seed has little bearing on the converged
#' objective in any case -- the greedy construction is a 2-approximation from
#' *any* start, and the drop-add tabu search determines the result. Ties break
#' to the smallest index, as in the matrix kernel, so when the two seed rules
#' coincide the whole trajectory matches the matrix path bit-for-bit. Pass
#' `DropAdd(seed=)` to override.
#'
#' @param colFn Column oracle; see [DropAdd()].
#' @param N Integer element count.
#' @param k Integer subset size (`2 <= k <= N`).
#' @param first Integer 1-based index of the seed element.
#' @return `.DropAddConstructColumn()` returns a `new.env(parent = emptyenv())`
#'   holding `S` (the selection, in add order), `inS`, `minDist`, `sumDist` and
#'   `minDistCount`, matching the record set of `.DropAddConstruct()`.
#' @references \insertAllCited{}
#' @keywords internal
.DropAddConstructColumn <- function(colFn, N, k, first) {
  st <- new.env(parent = emptyenv())
  st$S <- integer(k)
  st$S[1L] <- first
  st$inS <- logical(N)
  st$inS[first] <- TRUE

  col <- .DropAddColumn(colFn, first, N)
  st$minDist <- col
  st$sumDist <- col
  st$minDistCount <- rep(1L, N)
  st$minDist[first] <- Inf             # mask self
  st$minDistCount[first] <- 0L

  for (h in seq_len(k - 1L) + 1L) {
    xNew <- .DropAddPick(st)
    st$S[h] <- xNew
    st$inS[xNew] <- TRUE
    col <- .DropAddColumn(colFn, xNew, N)
    .DropAddApplyAdd(st, col, xNew)
    # xNew's own record, over the previously selected elements.
    dv <- col[st$S[seq_len(h - 1L)]]
    mn <- min(dv)
    st$minDist[xNew] <- mn
    st$minDistCount[xNew] <- sum(dv == mn)
  }
  # Return:
  st
}

#' DropAdd tabu search from a distance-column oracle
#'
#' The pure-R counterpart of `DropAdd_cpp()`, substituting an on-demand
#' `colFn(i)` call for the matrix-column read `dmat[, i]`. Implements the
#' distance-column oracle path of [DropAdd()] (dispatched there when `d` is a
#' function); see that function's *Distance-column oracle* section for the
#' user-facing contract. Every record update, tie-break and stopping rule
#' matches `src/dropadd.cpp` statement for statement, so on a symmetric oracle
#' the whole search trajectory -- not merely the final subset -- is identical to
#' the matrix path given the same seed.
#'
#' @inheritParams .DropAddConstructColumn
#' @param plateau Integer: stop after this many consecutive non-improving
#'   iterations.
#' @param maxSeconds Numeric wall-clock budget. Checked once per iteration (an
#'   iteration is dominated by two oracle calls, so `proc.time()` is free by
#'   comparison and the 1024-iteration stride of the C++ kernel is unnecessary);
#'   the check is skipped entirely when the budget is infinite, keeping the
#'   default path deterministic.
#' @param maxIter Integer cap on main-loop iterations.
#' @param trace Logical: also return the dropped/added index sequences, for the
#'   trajectory-identity tests.
#' @return `.DropAddFromColumn()` returns a list with the same shape as
#'   `DropAdd_cpp()`'s: `indices` (1-based, in FIFO buffer order), `objective`,
#'   `secondary`, `iters`, and -- when `trace` is `TRUE` -- `drops` and `adds`.
#' @keywords internal
.DropAddFromColumn <- function(colFn, N, k, first, plateau = 5000L,
                               maxSeconds = Inf,
                               maxIter = .Machine$integer.max,
                               trace = FALSE) {
  eps <- 1e-9
  st <- .DropAddConstructColumn(colFn, N, k, first)

  # -- Initial objective ---------------------------------------------------
  # `sum()` accumulates in long double, matching the C++'s `long double`
  # accumulator; halving is exact, so `secondary` is bit-identical either way.
  bestS       <- st$S
  bestMaxmin  <- min(st$minDist[st$S])
  bestSumpair <- sum(st$sumDist[st$S]) * 0.5
  bestScore   <- bestMaxmin + eps * bestSumpair

  head <- 1L                           # 1-based: FIFO drop position
  itersDone <- 0L
  noImprove <- 0L
  drops <- integer(0)
  adds  <- integer(0)
  timed <- is.finite(maxSeconds)
  t0 <- proc.time()[[3L]]

  # k == N leaves Add X(k) empty once x# is excluded, so no drop-add move
  # exists; the construction already selected every element.
  while (itersDone < maxIter && k < N) {
    if (noImprove >= plateau) {
      break
    }
    if (timed && proc.time()[[3L]] - t0 >= maxSeconds) {
      break
    }

    # 1. DROP: x# = S[head] (FIFO via circular buffer).
    xHash <- st$S[[head]]
    st$inS[xHash] <- FALSE
    survivors <- st$S[-head]           # S is rewritten only at the very end
    dXhash <- .DropAddColumn(colFn, xHash, N)
    st$sumDist <- st$sumDist - dXhash
    # Elements that had x# as a nearest peer lose one; those left with none need
    # their record rebuilt. (`<=`, not `==`, mirroring the kernel.)
    touched <- dXhash <= st$minDist
    touched[xHash] <- FALSE
    st$minDistCount[touched] <- st$minDistCount[touched] - 1L
    needRecompute <- which(touched & st$minDistCount == 0L)

    for (xx in needRecompute) {
      # Inverted relative to src/dropadd.cpp, which scans the k - 1 surviving
      # *columns* S[j]: from an oracle that would be k - 1 calls per triggering
      # iteration, the one place a literal port would be quadratic in oracle
      # calls. Reading xx's own column instead costs one call per entry of
      # needRecompute. min and the equality count are order-independent, and the
      # oracle is assumed symmetric, so the rebuilt records are bit-identical.
      peers <- survivors[survivors != xx]
      if (length(peers)) {
        dv <- .DropAddColumn(colFn, xx, N)[peers]
        mn <- min(dv)
        st$minDist[xx] <- mn
        st$minDistCount[xx] <- sum(dv == mn)
      } else {
        # Reachable only at k == 2, where x#'s sole peer is xx itself.
        st$minDist[xx] <- Inf
        st$minDistCount[xx] <- 0L
      }
    }

    # x#'s own record: distances to the surviving peers. With k >= 2 and finite
    # distances there is always at least one, so the kernel's unreachable
    # infinite-minimum branch has no counterpart here.
    dv <- dXhash[survivors]
    mn <- min(dv)
    st$minDist[xHash] <- mn
    st$minDistCount[xHash] <- sum(dv == mn)

    # 2. ADD over Add X(k) = Z - X(k), with x# barred for this iteration.
    xNew <- .DropAddPick(st, exclude = xHash)
    st$inS[xNew] <- TRUE
    colNew <- .DropAddColumn(colFn, xNew, N)
    .DropAddApplyAdd(st, colNew, xNew)
    dv <- colNew[survivors]
    mn <- min(dv)
    st$minDist[xNew] <- mn
    st$minDistCount[xNew] <- sum(dv == mn)

    # Write xNew into the head slot and advance the FIFO.
    st$S[head] <- xNew
    head <- if (head == k) 1L else head + 1L

    # 3. Test improvement of the best-known MMDPo solution.
    cm <- min(st$minDist[st$S])
    cs <- sum(st$sumDist[st$S]) * 0.5
    curScore <- cm + eps * cs
    if (curScore > bestScore) {
      bestS       <- st$S
      bestMaxmin  <- cm
      bestSumpair <- cs
      bestScore   <- curScore
      noImprove   <- 0L
    } else {
      noImprove <- noImprove + 1L
    }

    itersDone <- itersDone + 1L
    if (trace) {
      drops[[itersDone]] <- xHash
      adds[[itersDone]]  <- xNew
    }
  }

  out <- list(indices = bestS, objective = bestMaxmin,
              secondary = bestSumpair, iters = itersDone)
  if (trace) {
    out$drops <- drops
    out$adds  <- adds
  }
  # Return:
  out
}

# Internal test scaffolding (not part of the public API). Runs the production
# C++ DropAdd with tracing enabled and returns the dropped/added index
# sequences alongside the result, so the FIFO + tabu invariants can be asserted
# without exposing a `.trace` argument on `DropAdd()`. Kept deliberately thin:
# it mirrors only the `d`-path coercion and the single C++ call.
.DropAddTrace <- function(d, k, maxIter = NULL, plateau = 5000L,
                          maxSeconds = Inf) {
  dmat <- .AsDistMatrix(d)
  cppMaxIter <- if (is.null(maxIter)) .Machine$integer.max else as.integer(maxIter)
  out <- DropAdd_cpp(dmat, as.integer(k), as.double(maxSeconds),
                       cppMaxIter, as.integer(plateau), TRUE)
  list(
    indices = sort(as.integer(out$indices)),
    score   = as.numeric(out$objective),
    iters   = as.integer(out$iters),
    drops   = as.integer(out$drops),
    adds    = as.integer(out$adds)
  )
}

# ----- main entry point -----------------------------------------------------

#' DropAdd Tabu Search for the Max-Min Diversity Problem
#'
#' `DropAdd()` selects a maximally-dispersed subset of `k` points using the
#' DropAdd tabu search algorithm, which comprises a greedy construction
#' followed by a first-in, first-out drop-add tabu search, with streamlined
#' neighbour-evaluation tricks
#' \insertCite{@algorithms 1--4 in @Porumbel2011}{MaxMin}.
#'
#' @param k Integer; subset size, \eqn{2 \le k \le N}.
#' @param d A \code{dist} object, a square symmetric numeric matrix, or a
#'  distance-column function that takes an index `i` and returns the distances
#'  from `i` to each other element (optionally including the self-distance); see
#'  the *Distance-column oracle* section.
#' @param N Integer: the total number of elements. Required (and used) only on
#'  the distance-column oracle path, where it cannot be inferred from the
#'  closure; ignored for the matrix and coordinate paths.
#' @param points A numeric \eqn{N \times \mathrm{dim}} coordinate matrix (or an
#'  object coercible to one via \code{as.matrix}).
#'  Must be complete (no \code{NA}).
#'  Ignored if `d` specified.
#'  Avoids creating an \eqn{N \times N} distance matrix, enabling use at
#'  \eqn{N \ge 46340}).
#' @param plateau Integer; stop after this many consecutive drop-add
#'  iterations do not improve the score.
#' @param maxSeconds Numeric: terminate search after this many seconds have
#' elapsed.
#' @param seed Optional integer: a 1-based start index that overrides the
#'  construction's default warm-start seed (max-row-sum on the `d` matrix path,
#'  centroid-peripheral on the `points` path, two-sweep peripheral on the
#'  distance-column path), mirroring [FarFirst()]'s integer
#'  `strategy`. `NULL` (default) keeps the method's own seed. Not supported when
#'  candidate thinning binds (pass `maxCandidates = 0L`).
#' @templateVar default `46340L`
#' @templateVar default_basis the dense-distance-matrix feasibility ceiling
#'   (`floor(sqrt(.Machine$integer.max))`) the `points` path already crosses;
#'   below it nothing changes, at or above it the candidates are thinned to this
#'   size on the matrix-free path (no \eqn{m \times m} matrix is built)
#' @template maxCandidates
#' @templateVar progress_shows status messages are shown
#' @template progress
#'
#' @return `DropAdd()` returns an integer vector of length \code{k} containing the 1-based selected
#'   indices **sorted ascending** (unlike [FarFirst()], which returns
#'   farthest-first order), with attributes:
#'   \describe{
#'     \item{score}{numeric(1), achieved MaxMin objective
#'       \eqn{\min_{i \ne j \in S} d_{ij}}.}
#'     \item{secondary}{numeric(1), achieved sum of pairwise distances over
#'       \eqn{S} (upper-triangle sum).}
#'     \item{time_s}{numeric(1), wall-clock seconds spent.}
#'     \item{iters}{integer(1), main-loop iterations executed (excluding the
#'       construction phase).}
#'   }
#'   The vector has class `"MaxMinSelection"` and prints as a one-line summary
#'   (see [print.MaxMinSelection()]); it is otherwise an ordinary integer vector.
#'
#' @section Distance-column oracle:
#' When `d` is a **function**, it is treated as a distance-column oracle:
#' `d(i)` must return the distances from element `i` to every element (length
#' `N`, the self-distance ignored) or to every *other* element (length
#' `N - 1`, in index order). Supply `N`. This suits metrics with neither a
#' stored matrix nor a coordinate embedding -- on-demand tree-to-tree distances,
#' say -- and never materialises the \eqn{N \times N} matrix: memory is
#' \eqn{O(N)}.
#'
#' The search is the same algorithm against a different distance source, so on a
#' symmetric oracle the trajectory and result are identical to the matrix path
#' given the same seed. Two things differ:
#'
#' - **The seed.** The matrix path's max-row-sum warm start would need all `N`
#'   columns, i.e. the whole pairwise set. The oracle path substitutes a
#'   two-sweep peripheral seed costing two calls (see
#'   [.DropAddConstructColumn()]); `seed=` overrides it.
#' - **`maxCandidates` thinning is unavailable** (it builds a coreset distance
#'   matrix), and a binding cap warns rather than silently thinning.
#'
#' Cost is counted in **oracle calls**, not flops: `k` for the construction (two
#' more for the default seed), plus two per main-loop iteration (and one per
#' element whose nearest selected
#' peer just vanished). `plateau` and `maxSeconds` therefore set the distance
#' bill directly, and with the default `plateau` a run can request far more
#' columns than the `N` a full matrix would cost. Write the closure accordingly:
#' hoist everything reusable into the enclosing environment so a call pays only
#' for the `N` comparisons, and -- when `N` columns of cache are affordable --
#' memoise, which caps the oracle at `N` distinct evaluations.
#'
#' @references \insertAllCited{}
#'
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(200), ncol = 2)
#' DropAdd(5L, dist(pts))
#'
#' # Composable coreset: thin to 40 candidates with farthest-first, then run
#' # DropAdd on the coreset. Returned indices are original-space row indices.
#' suppressWarnings(DropAdd(5L, points = pts, maxCandidates = 40L))
#'
#' # Disable thinning (run on the full problem):
#' DropAdd(5L, points = pts, maxCandidates = 0L)
#'
#' # Distance-column oracle: `d` is a function of one index. Build it in a
#' # factory so everything the metric reuses is derived ONCE, in the enclosing
#' # environment, and each call only looks it up -- a closure that re-derives its
#' # representation per call pays that cost thousands of times over. `cache`
#' # memoises the columns, capping the oracle at one evaluation per element.
#' data("USArrests")
#' ArrestOracle <- function(dat) {
#'   scaled <- scale(as.matrix(dat))                 # derived once
#'   cache <- new.env(parent = emptyenv())
#'   function(i) {
#'     key <- as.character(i)
#'     if (is.null(cache[[key]])) {
#'       cache[[key]] <- sqrt(rowSums(sweep(scaled, 2, scaled[i, ], "-") ^ 2))
#'     }
#'     cache[[key]]
#'   }
#' }
#' arrests <- USArrests[, c("Murder", "Assault", "Rape")]
#' idx <- DropAdd(4L, ArrestOracle(arrests), N = nrow(arrests), plateau = 200L)
#' USArrests[idx, ]
#' @export
DropAdd <- function(k, d = NULL, plateau = 5000L, maxSeconds = Inf,
                    points = NULL, maxCandidates = 46340L, seed = NULL,
                    N = NULL) {
  progress <- getOption("MaxMin.progress", interactive())
  if (!is.null(points) && !is.null(d)) {
    stop("supply `d` or `points`, not both")
  }
  # Distance-column oracle path: `d` is a closure returning one matrix column at
  # a time, for metrics with neither a stored matrix nor a coordinate embedding
  # (e.g. on-demand tree-to-tree distances). Detected before .AsDistMatrix(),
  # which rejects a function, and before the thinning block, which would need a
  # coreset distance matrix. Mirrors FarFirst()'s dispatch.
  useOracle <- is.function(d)
  usePoints <- !useOracle && !is.null(points)
  if (useOracle) {
    if (is.null(N)) {
      stop("`N` (the element count) must be supplied when `d` is a ",
           "distance-column function")
    }
    N <- as.integer(N)
    if (length(N) != 1L || is.na(N) || N < 1L) {
      stop("`N` must be a single positive integer")
    }
    n <- N
  } else {
    if (!is.null(N)) {
      warning("`N` is ignored on the matrix/coordinate path; ",
              "it is only needed when `d` is a distance-column function")
    }
    if (usePoints) {
      points <- .AsPointsMatrix(points)
      n <- nrow(points)
    } else {
      dmat <- .AsDistMatrix(d)
      n <- nrow(dmat)
    }
  }
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L || k > n) {
    stop("`k` must be a single integer with 2 <= k <= n")
  }
  plateau <- as.integer(plateau)
  if (length(plateau) != 1L || is.na(plateau) ||
      plateau < 1L) {
    stop("`plateau` must be a single positive integer")
  }
  if (!is.numeric(maxSeconds) || length(maxSeconds) != 1L ||
      is.na(maxSeconds) || maxSeconds <= 0) {
    stop("`maxSeconds` must be a single positive numeric (or Inf)")
  }
  # `seed` overrides the construction's default warm-start (max-row-sum on the
  # matrix path, centroid-peripheral on the points path, two-sweep peripheral on
  # the distance-column path) with an explicit 1-based start index, mirroring
  # FarFirst()'s integer `strategy`. NULL keeps the default.
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.finite(seed) || seed < 1L || seed > n)) {
    stop("`seed` must be a single index in [1, n], or NULL for the default")
  }
  seed0 <- if (is.null(seed)) -1L else as.integer(seed) - 1L

  # Composable-coreset thinning: when `maxCandidates` binds, farthest-first
  # reduces the n candidates to an m-point coreset and DropAdd runs on that.
  # The recursive call passes `maxCandidates = 0L`, so the coreset is solved by
  # the real kernel exactly once; on the points source the subproblem stays on
  # the matrix-free path (no m x m matrix is built).
  mc <- .ResolveCap(maxCandidates, n, k)
  if (!is.na(mc) && useOracle) {
    # Thinning restricts the problem to an m-point coreset and hands the solver
    # the m x m submatrix, which a column oracle cannot supply without extra
    # machinery. Warn rather than silently substituting, as FarFirst() does for
    # its unreachable seeding strategies.
    warning("`maxCandidates` thinning is not supported on the distance-column ",
            "path; running on the full problem", call. = FALSE)
    mc <- NA_integer_
  }
  if (!is.na(mc)) {
    if (!is.null(seed)) {
      stop("`seed=` is not supported with candidate thinning; pass ",
           "`maxCandidates = 0L` to run on the full problem")
    }
    return(.FarFirstThin(
      k, mc,
      d      = if (usePoints) NULL else dmat,
      points = if (usePoints) points else NULL,
      RunOnSubset = function(d, points)
        DropAdd(k, d = d, points = points, plateau = plateau,
                maxSeconds = maxSeconds, maxCandidates = 0L),
      label = "DropAdd"))
  }

  t0 <- proc.time()[[3L]]
  eps <- 1e-9

  # --- Distance-column oracle path ------------------------------------------
  if (useOracle) {
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (column oracle; n = {n}, k = {k}, budget = {maxSeconds}s)",
        .auto_close = FALSE
      )
    }
    first <- if (is.null(seed)) .PeripheralSeedColumn(d, n) else as.integer(seed)
    out <- .DropAddFromColumn(d, n, k, first, plateau = plateau,
                              maxSeconds = maxSeconds)
    timeS <- proc.time()[[3L]] - t0
    if (progress) {
      itersMsg <- as.integer(out$iters)
      tkMsg    <- as.numeric(out$objective)
      cli::cli_process_done(
        msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
      )
    }
    return(.AsMaxMinSelection(structure(
      sort(as.integer(out$indices)),
      score     = as.numeric(out$objective),
      secondary = as.numeric(out$secondary),
      time_s    = timeS,
      iters     = as.integer(out$iters)
    ), "DropAdd"))
  }

  # --- Matrix-free (coordinate) path ----------------------------------------
  if (usePoints) {
    if (progress) {
      cli::cli_process_start(
        "DropAdd tabu search (n = {n}, k = {k}, budget = {maxSeconds}s)",
        .auto_close = FALSE
      )
    }
    out <- DropAdd_points_cpp(points, k, as.double(maxSeconds),
                                .Machine$integer.max, plateau, FALSE, seed0)
    timeS <- proc.time()[[3L]] - t0
    if (progress) {
      itersMsg <- as.integer(out$iters)
      tkMsg    <- as.numeric(out$objective)
      cli::cli_process_done(
        msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
      )
    }
    return(.AsMaxMinSelection(structure(
      sort(as.integer(out$indices)),
      score     = as.numeric(out$objective),
      secondary = as.numeric(out$secondary),
      time_s    = timeS,
      iters     = as.integer(out$iters)
    ), "DropAdd"))
  }

  # --- C++ drop-add tabu search (the sole compute path) -------------------
  if (progress) {
    cli::cli_process_start(
      "DropAdd tabu search (n = {n}, k = {k}, budget = {maxSeconds}s)",
      .auto_close = FALSE
    )
  }
  out <- DropAdd_cpp(dmat, k, as.double(maxSeconds),
                       .Machine$integer.max, plateau, FALSE, seed0)
  timeS <- proc.time()[[3L]] - t0
  if (progress) {
    itersMsg <- as.integer(out$iters)
    tkMsg    <- as.numeric(out$objective)
    cli::cli_process_done(
      msg = "DropAdd: {itersMsg} iters, T_k = {signif(tkMsg, 4)}, {round(timeS, 1)}s"
    )
  }
  .AsMaxMinSelection(structure(
    sort(as.integer(out$indices)),
    score     = as.numeric(out$objective),
    secondary = as.numeric(out$secondary),
    time_s    = timeS,
    iters     = as.integer(out$iters)
  ), "DropAdd")
}
