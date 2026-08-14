# seed.R
#
# Peripheral seeding strategies for Gonzalez farthest-first selection, and the
# ensemble driver that runs several and keeps the best by MinDist(). The single
# anchors are exposed through PickPoint(); FarFirst(strategy = ) selects among them
# (or the ensemble) and runs the greedy pass.

# Seeds available to the ensemble drivers. The matrix path lacks coordinates, so
# `"anti_centroid"` is reachable only from the coordinate path.
.kMatrixEnsembleSeeds <- c("peripheral", "random_furthest", "diameter",
                           "anti_medoid", "medoid", "rowsum", "rownorm")
.kPointEnsembleSeeds  <- c("anti_centroid", .kMatrixEnsembleSeeds)

#' Draw distinct furthest-point seeds from random pivots
#'
#' Used by [.GonzEnsemble()] and [.GonzEnsembleFromPoints()] to expand the
#' `"random_furthest"` token (see [FarFirst()]): walks distinct random pivots
#' (a partial shuffle of `1:nPts`, so no pivot is ever tried twice), resolves
#' each pivot's furthest-point seed via `SeedFunc`, and collects distinct seed
#' indices until `nSeeds` are found or the draw budget is spent. A `maxDraws`
#' cap bounds the work when the reachable seed pool is smaller than `nSeeds`.
#' Returns between 1 and `nSeeds` distinct indices (ascending); set a seed
#' (`set.seed()`) for a reproducible set.
#' @param SeedFunc Function mapping a pivot index to its furthest-point seed index.
#' @param nPts Integer number of points.
#' @param nSeeds Integer target number of distinct seeds (`>= 1`).
#' @param maxDraws Integer cap on the number of distinct pivots tried.
#'   Default `max(40 * nSeeds, 100)`.
#' @return `.DrawDistinctSeeds()` returns an integer vector of distinct seed indices (length in `[1, nSeeds]`).
#' @keywords internal
.DrawDistinctSeeds <- function(SeedFunc, nPts, nSeeds, maxDraws = NULL) {
  nSeeds <- as.integer(nSeeds)
  if (is.null(maxDraws)) maxDraws <- max(40L * nSeeds, 100L)
  seen  <- logical(nPts)
  found <- 0L
  for (r in sample.int(nPts, min(as.integer(maxDraws), nPts))) {
    s <- SeedFunc(r)
    if (!seen[s]) {
      seen[s] <- TRUE
      found <- found + 1L
      if (found >= nSeeds) break
    }
  }
  # Return:
  which(seen)
}

#' Squared distance of every point to the coordinate anti_centroid
#'
#' The `O(N * dim)` basis of the `"anti_centroid"` seed: its argmax is the point
#' farthest from the coordinate mean, an approximate diameter endpoint.
#' @param points A `double` `N x dim` coordinate matrix.
#' @return `.CentroidSqDist()` returns a numeric vector of length `N` of squared distances to the mean.
#' @keywords internal
.CentroidSqDist <- function(points) {
  mu  <- colMeans(points)
  dev <- points - rep(mu, each = nrow(points))
  rowSums(dev * dev)
}


#' Expand ensemble anchor names into labelled seed specs
#'
#' Maps each anchor name to a `list(label, s1)` spec. The `"random_furthest"`
#' token expands to one spec per element of `rfSeeds` (already-resolved seed
#' indices, labelled `random_furthest1`, ...); an empty vector contributes none.
#' @param anchors Character vector of (de-duplicated) anchor names.
#' @param rfSeeds Integer vector of already-resolved furthest-point seed indices
#'   for the `"random_furthest"` token.
#' @param anchorSeed Function mapping a deterministic anchor name to an integer
#'   seed index.
#' @return `.ExpandAnchors()` returns a list of `list(label, s1)` specs.
#' @keywords internal
.ExpandAnchors <- function(anchors, rfSeeds, anchorSeed) {
  specs <- list()
  for (name in anchors) {
    if (name == "random_furthest") {
      for (j in seq_along(rfSeeds)) {
        specs[[length(specs) + 1L]] <- list(
          label = paste0("random_furthest", j),
          s1    = rfSeeds[[j]]
        )
      }
    } else {
      sd <- anchorSeed(name)
      specs[[length(specs) + 1L]] <- list(label = name, s1 = sd)
    }
  }
  if (length(specs) == 0L) {
    stop("no seed strategies to run")
  }
  specs
}

#' Resolve an expanded ensemble into the winning subset
#'
#' Shared tail of the two ensemble drivers: solves the specs' distinct seeds in
#' one batch via the driver's `RunPasses` closure, then returns the subset
#' maximising \eqn{T_k}. The returned vector carries the `strategy_results` (one
#' record per label) and `winning_strategy` (all tied-best labels) attributes.
#' @param expanded List of `list(label, s1)` specs from [.ExpandAnchors()].
#' @param labels Character vector of labels (one per spec).
#' @param RunPasses Closure mapping a vector of distinct seeds to a list of
#'   `list(idx, tK)`, one per seed and in the same order.
#' @return `.ResolveEnsemble()` returns an integer vector of selected indices with attributes.
#' @keywords internal
.ResolveEnsemble <- function(expanded, labels, RunPasses) {
  # Two anchors can resolve to the same seed (and a random draw can land on an
  # anchor's), so the batch is de-duplicated before dispatch and mapped back
  # after: each distinct seed is solved exactly once, and every label still
  # gets its own record.
  seeds    <- vapply(expanded, function(e) as.integer(e$s1), integer(1L))
  distinct <- unique(seeds)
  runs     <- RunPasses(distinct)
  back     <- match(seeds, distinct)
  strategyResults <- lapply(seq_along(expanded), function(i) {
    g <- runs[[back[[i]]]]
    list(s1 = expanded[[i]]$s1, idx = g$idx, t_k = g$tK)
  })
  names(strategyResults) <- labels

  # Maximise T_k (the selection's min pairwise distance; larger = better spread).
  # which.max() skips NA and first-wins on ties, matching the old hand-rolled
  # loop; all-NA (e.g. n == 1) falls back to the first strategy.
  tks     <- vapply(strategyResults, `[[`, numeric(1L), "t_k")
  allNa   <- all(is.na(tks))
  bestI   <- if (allNa) 1L else which.max(tks)
  bestTk  <- tks[[bestI]]
  # When every T_k is NA (e.g. n == 1, each pass selects a single point) all
  # strategies are trivially tied-best, so report them all -- matching the
  # documented "all tied-best strategies" contract (F-602). Otherwise the
  # winners are the strategies achieving the best finite T_k.
  winners <- if (allNa) seq_along(tks) else which(tks == bestTk)

  result <- strategyResults[[bestI]]$idx
  # Expose the winning T_k as `score`, matching the bare FarFirst() pass and the
  # DropAdd()/Grasp() return contract.
  attr(result, "score")            <- bestTk
  attr(result, "strategy_results") <- strategyResults
  attr(result, "winning_strategy") <- labels[winners]
  result
}

#' Peripheral seed index for Gonzalez selection (distance matrix)
#'
#' @param d Square numeric distance matrix.
#' @param strategy Anchor name; see [PickPoint()]. Also accepts `"first"` (1).
#' @return `.PickPoint()` returns an integer seed index.
#' @keywords internal
.PickPoint <- function(d, strategy) {
  switch(strategy,
    first   = 1L,
    anti_centroid = stop("`anti_centroid` strategy requires coordinates; supply `points=` ",
                    "or use `peripheral` on the distance-matrix path"),
    medoid  = as.integer(which.min(rowSums(d))),
    rowsum  = as.integer(which.max(rowSums(d))),
    # Fused C++ scan: rowSums(d ^ 2) materialises the full squared matrix.
    rownorm = as.integer(which.max(RowSqSumsFromMatrix_cpp(d, .NThreads()))),
    anti_medoid = {
      med <- which.min(rowSums(d))
      dd  <- d[, med]
      dd[med] <- -Inf
      as.integer(which.max(dd))
    },
    diameter = {
      # C++ off-diagonal first-max scan: the R idiom copies the full matrix
      # to plant -Inf on the diagonal, then scans it twice more.
      dm <- MatrixOffDiagMax_cpp(d, .NThreads())
      if (!is.finite(dm[[1L]]) || dm[[1L]] <= 0) {
        1L
      } else {
        as.integer(dm[[2L]])
      }
    },
    peripheral = {
      s1 <- which.max(d[, 1L])
      as.integer(which.max(d[, s1]))
    },
    random_furthest = {
      r <- sample.int(nrow(d), 1L)
      as.integer(which.max(d[, r]))
    },
    stop("Unknown seed strategy: ", strategy)
  )
}

#' Peripheral seed index for Gonzalez selection (coordinates)
#'
#' Coordinate counterpart of [.PickPoint()]; each anchor is computed from the
#' `…FromPoints_cpp` primitives, bit-identical to the matrix path on Euclidean
#' data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param strategy Anchor name; see [PickPoint()]. Also accepts `"first"` (1).
#' @return `.PickPoints()` returns an integer seed index.
#' @keywords internal
.PickPoints <- function(points, strategy) {
  switch(strategy,
    first   = 1L,
    anti_centroid = as.integer(which.max(.CentroidSqDist(points))),
    medoid  = as.integer(which.min(RowSumsFromPoints_cpp(points, .NThreads()))),
    rowsum  = as.integer(which.max(RowSumsFromPoints_cpp(points, .NThreads()))),
    rownorm = as.integer(which.max(RowSqSumsFromPoints_cpp(points, .NThreads()))),
    anti_medoid = {
      med <- which.min(RowSumsFromPoints_cpp(points, .NThreads()))
      dd  <- EuclidColFromPoints_cpp(points, med)
      dd[med] <- -Inf
      as.integer(which.max(dd))
    },
    diameter = {
      diam <- DiameterFromPoints_cpp(points, .NThreads())
      if (!is.finite(diam[1L]) || diam[1L] <= 0) {
        1L
      } else {
        as.integer(diam[2L])
      }
    },
    peripheral = {
      s1 <- which.max(EuclidColFromPoints_cpp(points, 1L))
      as.integer(which.max(EuclidColFromPoints_cpp(points, s1)))
    },
    random_furthest = {
      r <- sample.int(nrow(points), 1L)
      as.integer(which.max(EuclidColFromPoints_cpp(points, r)))
    },
    stop("Unknown seed strategy: ", strategy)
  )
}

#' Seed to initialize farthest-first selection
#'
#' `PickPoint()` implements a range of strategies to select a seed for greedy
#' farthest-first selection. Propitious seeds yield better solutions.
#'
#'
#' @param d A `dist` object or square symmetric numeric matrix. Ignored when
#'   `points` is supplied.
#' @param points Optional `N x dim` numeric coordinate matrix; when supplied the
#'   seed is computed from coordinates in `O(N)` memory. Required for the
#'   `"anti_centroid"` anchor, which has no distance-matrix form.
#' @param strategy Character specifying method to employ:
#'
#' \describe{
#'   \item{`"peripheral"`(default)}{Two sweeps: the point furthest from point 1,
#'    then the point furthest from that (a diameter-endpoint approximation).
#'    \eqn{O(N)}.}
#'   \item{`"anti_centroid"`}{The point farthest from the coordinate mean
#'     (\eqn{\argmax ||x - x_bar||}). \eqn{O(N * dim)}. Requires `points`.}
#'   \item{`"random_furthest"`}{The point furthest from a random pivot.
#'   \eqn{O(N)}.}
#'   \item{`"diameter"`}{A row endpoint of the diameter pair (the maximum
#'     pairwise distance).}
#'   \item{`"medoid"`}{The 1-median (medoid): the point minimising the sum of
#'     distances to all others.}
#'   \item{`"anti_medoid"`}{The point furthest from the 1-median (medoid).}
#'   \item{`"rowsum"`}{The point maximising the sum of distances to all others
#'     (the 1-anti-median).}
#'   \item{`"rownorm"`}{The point maximising \eqn{\sqrt(\sum{d^2})}, the L2
#'    counterpart of `"rowsum"`.}
#' }
#'
#' @return `PickPoint()` returns an integer that identifies the index of a
#' proposed seed in `d` or `points`.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' PickPoint(d, strategy = "diameter")
#' FarFirst(5L, d, strategy = PickPoint(d, strategy = "diameter"))
#' @seealso [FarFirst()], which seeds and runs the greedy pass in one call.
#' @export
PickPoint <- function(d = NULL, points = NULL,
                       strategy = c("peripheral", "anti_centroid",
                                    "random_furthest", "diameter",
                                    "anti_medoid", "medoid", "rowsum",
                                    "rownorm")) {
  strategy <- match.arg(strategy)
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    return(.PickPoints(points, strategy))
  }
  d <- .AsDistMatrix(d)
  .PickPoint(d, strategy)
}

#' Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)
#'
#' Runs Gonzalez from each requested peripheral anchor and returns the subset
#' maximising \eqn{T_k}. Internal driver for the ensemble path of [FarFirst()]
#' (triggered when `strategy` is a character vector of length > 1 or
#' `"random_furthest"`). The `"random_furthest"` token draws `nSeeds` distinct
#' furthest-point seeds via [.DrawDistinctSeeds()]. The returned vector carries
#' `strategy_results` and `winning_strategy` (character vector of all tied-best
#' strategies, with random starts labelled `random_furthest1`,
#' `random_furthest2`, ...) attributes.
#' @param d Square numeric distance matrix (already coerced).
#' @param m Integer subset size (`1 <= m < nrow(d)`).
#' @param anchors Character vector of anchor names.
#' @param nSeeds Integer number of distinct random-furthest seeds to draw when
#'   `"random_furthest"` is in `anchors`.
#' @return `.GonzEnsemble()` returns an integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsemble <- function(d, m, anchors = "peripheral", nSeeds = 3L) {
  d <- .AsDistMatrix(d)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 0L) {
    stop("`m` must be a single non-negative integer")
  }
  if (is.null(anchors) || length(anchors) == 0L) {
    stop("`anchors` must name at least one strategy")
  }
  anchors <- unique(match.arg(
    anchors,
    choices = .kMatrixEnsembleSeeds,
    several.ok = TRUE
  ))
  nPts <- nrow(d)
  if (m >= nPts) return(structure(seq_len(nPts), score = NA_real_))
  if (m == 0L)   return(integer(0))

  lazy <- new.env(parent = emptyenv())
  GetRowSums <- function() {
    lazy$rowSums <- lazy$rowSums %||% rowSums(d)
    lazy$rowSums
  }
  GetRowSqSums <- function() {
    # Fused C++ scan: rowSums(d ^ 2) materialises the full squared matrix.
    lazy$rowSqSums <- lazy$rowSqSums %||% RowSqSumsFromMatrix_cpp(d, .NThreads())
    lazy$rowSqSums
  }
  GetMedoid <- function() {
    lazy$medoid <- lazy$medoid %||% as.integer(which.min(GetRowSums()))
    lazy$medoid
  }

  # One batched call: the kernel runs the passes one per thread and reports
  # each one's T_k (computed for free during its greedy pass), so nothing here
  # re-scores with a d[idx, idx] subset.
  RunPasses <- function(seeds) .MaximinMulti(m, seeds, d = d)

  AnchorSeed <- function(name) {
    switch(name,
      diameter = {
        # C++ off-diagonal first-max scan; see .PickPoint's diameter branch.
        dm <- MatrixOffDiagMax_cpp(d, .NThreads())
        if (!is.finite(dm[[1L]]) || dm[[1L]] <= 0) {
          1L
        } else {
          as.integer(dm[[2L]])
        }
      },
      anti_medoid = {
        med <- GetMedoid()
        distFromMedoid <- d[, med]
        distFromMedoid[med] <- -Inf
        as.integer(which.max(distFromMedoid))
      },
      medoid  = GetMedoid(),
      peripheral = {
        s1 <- which.max(d[, 1L])
        as.integer(which.max(d[, s1]))
      },
      rowsum  = as.integer(which.max(GetRowSums())),
      rownorm = as.integer(which.max(GetRowSqSums()))
    )
  }

  rfSeeds <- if ("random_furthest" %in% anchors) {
    .DrawDistinctSeeds(function(r) which.max(d[, r]), nPts, nSeeds)
  } else {
    integer(0)
  }
  expanded <- .ExpandAnchors(anchors, rfSeeds, AnchorSeed)
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  .ResolveEnsemble(expanded, labels, RunPasses)
}

#' Coordinate (matrix-free) multi-anchor Gonzalez ensemble
#'
#' Coordinate counterpart of [.GonzEnsemble()]; each anchor seed and the greedy
#' expansion are computed from `points` via the coordinate primitives, so the
#' returned indices and attributes match the matrix path on Euclidean data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param m Integer subset size.
#' @inheritParams .GonzEnsemble
#' @return `.GonzEnsembleFromPoints()` returns an integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsembleFromPoints <- function(points, m, anchors = "random_furthest",
                                    nSeeds = 3L) {
  points <- .AsPointsMatrix(points)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 0L) {
    stop("`m` must be a single non-negative integer")
  }
  if (is.null(anchors) || length(anchors) == 0L) {
    stop("`anchors` must name at least one strategy")
  }
  anchors <- unique(match.arg(
    anchors,
    choices = .kPointEnsembleSeeds,
    several.ok = TRUE
  ))
  nPts <- nrow(points)
  if (m >= nPts) return(structure(seq_len(nPts), score = NA_real_))
  if (m == 0L)   return(integer(0))

  lazy <- new.env(parent = emptyenv())
  # When the requested anchors span BOTH row-aggregate families, one fused
  # pair sweep fills both: each accumulator receives the identical summands
  # in the identical order as its dedicated kernel (bit-identical results),
  # while each pair's distance -- and its sqrt -- is computed once, not twice.
  if (any(c("anti_medoid", "medoid", "rowsum") %in% anchors) &&
        "rownorm" %in% anchors) {
    both <- RowSumsSqFromPoints_cpp(points, .NThreads())
    lazy$rowSums   <- both[[1L]]
    lazy$rowSqSums <- both[[2L]]
  }
  GetRowSums <- function() {
    lazy$rowSums <- lazy$rowSums %||% RowSumsFromPoints_cpp(points, .NThreads())
    lazy$rowSums
  }
  GetRowSqSums <- function() {
    lazy$rowSqSums <- lazy$rowSqSums %||% RowSqSumsFromPoints_cpp(points, .NThreads())
    lazy$rowSqSums
  }
  GetDiameter <- function() {
    lazy$diameter <- lazy$diameter %||% DiameterFromPoints_cpp(points, .NThreads())
    lazy$diameter
  }
  GetCentroidD2 <- function() {
    lazy$centroidD2 <- lazy$centroidD2 %||% .CentroidSqDist(points)
    lazy$centroidD2
  }
  GetMedoid <- function() {
    lazy$medoid <- lazy$medoid %||% as.integer(which.min(GetRowSums()))
    lazy$medoid
  }

  # One batched call: the kernel runs the passes one per thread and reports
  # each one's T_k (computed for free during its greedy pass), so nothing here
  # re-runs stats::dist() on the selected sub-coordinates.
  RunPasses <- function(seeds) .MaximinMulti(m, seeds, points = points)

  AnchorSeed <- function(name) {
    switch(name,
      diameter = {
        diam <- GetDiameter()
        if (!is.finite(diam[1L]) || diam[1L] <= 0) {
          1L
        } else {
          as.integer(diam[2L])
        }
      },
      anti_medoid = {
        med <- GetMedoid()
        distFromMedoid <- EuclidColFromPoints_cpp(points, med)
        distFromMedoid[med] <- -Inf
        as.integer(which.max(distFromMedoid))
      },
      anti_centroid  = as.integer(which.max(GetCentroidD2())),
      medoid    = GetMedoid(),
      peripheral = {
        s1 <- which.max(EuclidColFromPoints_cpp(points, 1L))
        as.integer(which.max(EuclidColFromPoints_cpp(points, s1)))
      },
      rowsum  = as.integer(which.max(GetRowSums())),
      rownorm = as.integer(which.max(GetRowSqSums()))
    )
  }

  rfSeeds <- if ("random_furthest" %in% anchors) {
    .DrawDistinctSeeds(
      function(r) which.max(EuclidColFromPoints_cpp(points, r)), nPts, nSeeds
    )
  } else {
    integer(0)
  }
  expanded <- .ExpandAnchors(anchors, rfSeeds, AnchorSeed)
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  .ResolveEnsemble(expanded, labels, RunPasses)
}
