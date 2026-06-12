# seed.R
#
# Peripheral seeding strategies for Gonzalez farthest-first selection, and the
# ensemble driver that runs several and keeps the best by MinDist(). The single
# anchors are exposed through MaxMinSeed(); FarFirst(strategy = ) selects among them
# (or the ensemble) and runs the greedy pass.

# The default seed ensemble: the `"random_furthest"` token alone, which draws
# `nseeds` (default 8) distinct furthest-point seeds via the session RNG and
# returns the best Gonzalez pass. Set a seed (`set.seed()`) for a reproducible
# draw. The deterministic anchors (`"centroid"`, `"peripheral"`, ...) remain
# available as opt-in `strategy=` options.
.kDefaultEnsemble <- "random_furthest"

# Seeds available to the ensemble drivers. The matrix path lacks coordinates, so
# `"centroid"` is reachable only from the coordinate path.
.kMatrixEnsembleSeeds <- c("peripheral", "random_furthest", "diameter",
                           "anti_medoid", "medoid", "rowsum", "rownorm")
.kPointEnsembleSeeds  <- c("centroid", .kMatrixEnsembleSeeds)

# Default number of distinct random-furthest seeds.
.kDefaultNSeeds <- 8L

#' Draw distinct furthest-point seeds from random pivots
#'
#' Used by [.GonzEnsemble()] and [.GonzEnsembleFromPoints()] to expand the
#' `"random_furthest"` token (see [FarFirst()]):
#' repeatedly draws a random pivot with the session RNG, resolves its
#' furthest-point seed via `seedFn`, and collects distinct seed indices until
#' `nseeds` are found. Two bounds stop the loop when the reachable seed pool is
#' smaller than `nseeds`: a consecutive-miss limit (the pool is likely exhausted
#' once many draws in a row yield only already-seen seeds) and an absolute draw
#' budget. Returns between 1 and `nseeds` distinct indices; set a seed
#' (`set.seed()`) for a reproducible set.
#' @param seedFn Function mapping a pivot index to its furthest-point seed index.
#' @param nPts Integer number of points.
#' @param nseeds Integer target number of distinct seeds (`>= 1`).
#' @param maxDraws Integer absolute draw budget. Default `max(40 * nseeds, 100)`.
#' @param missLimit Integer consecutive-miss limit. Default `max(8 * nseeds, 30)`.
#' @return Integer vector of distinct seed indices (length in `[1, nseeds]`).
#' @keywords internal
.DrawDistinctSeeds <- function(seedFn, nPts, nseeds, maxDraws = NULL,
                               missLimit = NULL) {
  nseeds <- as.integer(nseeds)
  if (is.null(maxDraws))  maxDraws  <- max(40L * nseeds, 100L)
  if (is.null(missLimit)) missLimit <- max(8L * nseeds, 30L)
  seen  <- integer(0)
  draws <- 0L
  miss  <- 0L
  while (length(seen) < nseeds && draws < maxDraws && miss < missLimit) {
    s <- as.integer(seedFn(sample.int(nPts, 1L)))
    draws <- draws + 1L
    if (s %in% seen) {
      miss <- miss + 1L
    } else {
      seen <- c(seen, s)
      miss <- 0L
    }
  }
  seen
}

#' Squared distance of every point to the coordinate centroid
#'
#' The `O(N * dim)` basis of the `"centroid"` seed: its argmax is the point
#' farthest from the coordinate mean, an approximate diameter endpoint.
#' @param points A `double` `N x dim` coordinate matrix.
#' @return Numeric vector of length `N` of squared distances to the mean.
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
#' @return List of `list(label, s1)` specs.
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
#' Shared tail of the two ensemble drivers: solves each expanded spec via the
#' driver's `RunGonz` closure (which de-duplicates repeated seeds through its
#' own cache), then returns the subset maximising \eqn{T_k}. The returned vector
#' carries the `strategy_results` (one record per label) and `winning_strategy`
#' (all tied-best labels) attributes.
#' @param expanded List of `list(label, s1)` specs from [.ExpandAnchors()].
#' @param labels Character vector of labels (one per spec).
#' @param RunGonz Closure mapping a seed `s1` to `list(idx, tK)`.
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.ResolveEnsemble <- function(expanded, labels, RunGonz) {
  strategyResults <- lapply(expanded, function(e) {
    g <- RunGonz(e$s1)
    list(s1 = e$s1, idx = g$idx, t_k = g$tK)
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
#' @param strategy Anchor name; see [MaxMinSeed()]. Also accepts `"first"` (1).
#' @return Integer seed index.
#' @keywords internal
.MaxMinSeed <- function(d, strategy) {
  switch(strategy,
    first   = 1L,
    centroid = stop("`centroid` strategy requires coordinates; supply `points=` ",
                    "or use `peripheral` on the distance-matrix path"),
    medoid  = as.integer(which.min(rowSums(d))),
    rowsum  = as.integer(which.max(rowSums(d))),
    rownorm = as.integer(which.max(rowSums(d ^ 2))),
    anti_medoid = {
      med <- which.min(rowSums(d))
      dd  <- d[, med]
      dd[med] <- -Inf
      as.integer(which.max(dd))
    },
    diameter = {
      dOff <- d
      diag(dOff) <- -Inf
      dMax <- max(dOff)
      if (!is.finite(dMax) || dMax <= 0) {
        1L
      } else {
        as.integer(arrayInd(which.max(dOff), dim(dOff))[1L, 1L])
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
#' Coordinate counterpart of [.MaxMinSeed()]; each anchor is computed from the
#' `…FromPoints_cpp` primitives, bit-identical to the matrix path on Euclidean
#' data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param strategy Anchor name; see [MaxMinSeed()]. Also accepts `"first"` (1).
#' @return Integer seed index.
#' @keywords internal
.MaxMinSeedPoints <- function(points, strategy) {
  switch(strategy,
    first   = 1L,
    centroid = as.integer(which.max(.CentroidSqDist(points))),
    medoid  = as.integer(which.min(RowSumsFromPoints_cpp(points))),
    rowsum  = as.integer(which.max(RowSumsFromPoints_cpp(points))),
    rownorm = as.integer(which.max(RowSqSumsFromPoints_cpp(points))),
    anti_medoid = {
      med <- which.min(RowSumsFromPoints_cpp(points))
      dd  <- EuclidColFromPoints_cpp(points, med)
      dd[med] <- -Inf
      as.integer(which.max(dd))
    },
    diameter = {
      diam <- DiameterFromPoints_cpp(points)
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

#' Peripheral seed index for Gonzalez farthest-first selection
#'
#' Returns the index of a single element, used to seed [FarFirst()].
#'
#' Anchors:
#' \describe{
#'   \item{`"centroid"`}{The point farthest from the coordinate mean
#'     (\eqn{argmax ||x - x_bar||}), an \eqn{O(N * dim)} approximate diameter
#'     endpoint.
#'     Computed from coordinates, so it requires `points`; it is unavailable on
#'     the distance-matrix path, where `"peripheral"` serves the same role.}
#'   \item{`"peripheral"`}{Two sweeps: the point furthest from point 1, then
#'     the point furthest from that (a diameter-endpoint approximation), in
#'     \eqn{O(N)}. The only anchor reachable from a distance-column oracle (the
#'     function path of [FarFirst()]).}
#'   \item{`"random_furthest"`}{The point furthest from a random pivot, in
#'     \eqn{O(N)}.}
#'   \item{`"diameter"`}{A row endpoint of the diameter pair (the maximum
#'     pairwise distance). Degenerate data (`d_max <= 0`) falls back to 1.}
#'   \item{`"anti_medoid"`}{The point furthest from the 1-median (medoid).}
#'   \item{`"medoid"`}{The 1-median.}
#'   \item{`"rowsum"`}{The point maximising the sum of distances to all others
#'     (the 1-anti-median).}
#'   \item{`"rownorm"`}{The point maximising \eqn{\sqrt(\sum{d^2})}, the L2
#'    counterpart of `"rowsum"`.}
#' }
#'
#' @param d A `dist` object or square symmetric numeric matrix. Ignored when
#'   `points` is supplied.
#' @param points Optional `N x dim` numeric coordinate matrix; when supplied the
#'   seed is computed from coordinates in `O(N)` memory. Required for the
#'   `"centroid"` anchor, which has no distance-matrix form.
#' @param strategy Anchor name (see above). Default `"peripheral"`.
#' @return Integer seed index in `[1, N]`.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' MaxMinSeed(d, strategy = "diameter")
#' FarFirst(d, 5L, strategy = MaxMinSeed(d, strategy = "diameter"))
#' @seealso [FarFirst()], which seeds and runs the greedy pass in one call.
#' @export
MaxMinSeed <- function(d = NULL, points = NULL,
                       strategy = c("peripheral", "centroid", "random_furthest",
                                    "diameter", "anti_medoid", "medoid", "rowsum",
                                    "rownorm")) {
  strategy <- match.arg(strategy)
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    return(.MaxMinSeedPoints(points, strategy))
  }
  d <- .AsDistMatrix(d)
  .MaxMinSeed(d, strategy)
}

#' Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)
#'
#' Runs Gonzalez from each requested peripheral anchor and returns the subset
#' maximising \eqn{T_k}. Internal driver for the ensemble path of [FarFirst()]
#' (triggered when `strategy` is a character vector of length > 1 or
#' `"random_furthest"`). The `"random_furthest"` token draws `nseeds` distinct
#' furthest-point seeds via [.DrawDistinctSeeds()]. The returned vector carries
#' `strategy_results` and `winning_strategy` (character vector of all tied-best
#' strategies, with random starts labelled `random_furthest1`,
#' `random_furthest2`, ...) attributes.
#' @param d Square numeric distance matrix (already coerced).
#' @param m Integer subset size (`1 <= m < nrow(d)`).
#' @param anchors Character vector of anchor names.
#' @param nseeds Integer number of distinct random-furthest seeds to draw when
#'   `"random_furthest"` is in `anchors`.
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsemble <- function(d, m, anchors = "peripheral",
                          nseeds = .kDefaultNSeeds) {
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
    lazy$rowSqSums <- lazy$rowSqSums %||% rowSums(d ^ 2)
    lazy$rowSqSums
  }
  GetDOffdiag <- function() {
    lazy$dOffdiag <- lazy$dOffdiag %||% { tmp <- d; diag(tmp) <- -Inf; tmp }
    lazy$dOffdiag
  }
  GetMedoid <- function() {
    lazy$medoid <- lazy$medoid %||% as.integer(which.min(GetRowSums()))
    lazy$medoid
  }

  gonzCache <- new.env(parent = emptyenv())
  RunGonz <- function(s1) {
    key <- as.character(s1)
    gonzCache[[key]] %||% {
      idx <- .MaximinFrom(d, m, first = s1)
      # The kernel computes T_k (min pairwise distance) for free during the
      # greedy pass; read it rather than re-scoring with a d[idx, idx] subset.
      tK  <- attr(idx, "t_k") %||% NA_real_
      attr(idx, "t_k") <- NULL
      res  <- list(idx = idx, tK = tK)
      gonzCache[[key]] <- res
      res
    }
  }

  AnchorSeed <- function(name) {
    switch(name,
      diameter = {
        dOff <- GetDOffdiag()
        dMax <- max(dOff)
        if (!is.finite(dMax) || dMax <= 0) {
          1L
        } else {
          as.integer(arrayInd(which.max(dOff), dim(dOff))[1L, 1L])
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
    .DrawDistinctSeeds(function(r) which.max(d[, r]), nPts, nseeds)
  } else {
    integer(0)
  }
  expanded <- .ExpandAnchors(anchors, rfSeeds, AnchorSeed)
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  .ResolveEnsemble(expanded, labels, RunGonz)
}

#' Coordinate (matrix-free) multi-anchor Gonzalez ensemble
#'
#' Coordinate counterpart of [.GonzEnsemble()]; each anchor seed and the greedy
#' expansion are computed from `points` via the coordinate primitives, so the
#' returned indices and attributes match the matrix path on Euclidean data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param m Integer subset size.
#' @inheritParams .GonzEnsemble
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsembleFromPoints <- function(points, m, anchors = .kDefaultEnsemble,
                                    nseeds = .kDefaultNSeeds) {
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
  GetRowSums <- function() {
    lazy$rowSums <- lazy$rowSums %||% RowSumsFromPoints_cpp(points)
    lazy$rowSums
  }
  GetRowSqSums <- function() {
    lazy$rowSqSums <- lazy$rowSqSums %||% RowSqSumsFromPoints_cpp(points)
    lazy$rowSqSums
  }
  GetDiameter <- function() {
    lazy$diameter <- lazy$diameter %||% DiameterFromPoints_cpp(points)
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

  gonzCache <- new.env(parent = emptyenv())
  RunGonz <- function(s1) {
    key <- as.character(s1)
    gonzCache[[key]] %||% {
      idx <- .MaximinFromPoints(points, m, first = s1)
      # The kernel computes T_k (min pairwise distance) for free during the
      # greedy pass; read it rather than re-running stats::dist() on the
      # selected sub-coordinates.
      tK  <- attr(idx, "t_k") %||% NA_real_
      attr(idx, "t_k") <- NULL
      res <- list(idx = idx, tK = tK)
      gonzCache[[key]] <- res
      res
    }
  }

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
      centroid  = as.integer(which.max(GetCentroidD2())),
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
      function(r) which.max(EuclidColFromPoints_cpp(points, r)), nPts, nseeds
    )
  } else {
    integer(0)
  }
  expanded <- .ExpandAnchors(anchors, rfSeeds, AnchorSeed)
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  .ResolveEnsemble(expanded, labels, RunGonz)
}
