# seed.R
#
# Peripheral seeding strategies for Gonzalez farthest-first selection, and the
# ensemble driver that runs several and keeps the best by MinDist(). The single
# anchors are exposed through MaxMinSeed(); FarFirst(seed = ) selects among them
# (or the ensemble) and runs the greedy pass.

# The default seed ensemble: the `"random_furthest"` token alone, which expands
# to one start per `pivots` element (three pivots drawn with the session RNG by
# default) -- a best-of-three cheap O(N) selection. Set a seed (`set.seed()`)
# for a reproducible draw. The deterministic anchors (`"centroid"`,
# `"peripheral"`, ...) remain available as opt-in `seed=` strategies.
.kDefaultEnsemble <- "random_furthest"

# Seeds available to the ensemble drivers. The matrix path lacks coordinates, so
# `"centroid"` is reachable only from the coordinate path.
.kMatrixEnsembleSeeds <- c("peripheral", "random_furthest", "diameter",
                           "anti_medoid", "medoid", "rowsum", "rownorm")
.kPointEnsembleSeeds  <- c("centroid", .kMatrixEnsembleSeeds)

# Number of random-furthest pivots drawn when `pivots` is left unspecified.
.kDefaultRandomStarts <- 3L

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
#' token expands to one spec per element of `pivots`, each seeded at the point
#' furthest from that pivot (labelled `random_furthest1`, ...); an empty
#' `pivots` contributes none.
#' @param anchors Character vector of (de-duplicated) anchor names.
#' @param pivots Integer vector of pivot indices the `"random_furthest"` token
#'   expands over (one start per pivot).
#' @param anchorSeed Function mapping a deterministic anchor name to an integer
#'   seed index.
#' @param rfSeed Function mapping a pivot index to the furthest-point seed.
#' @return List of `list(label, s1)` specs.
#' @keywords internal
.ExpandAnchors <- function(anchors, pivots, anchorSeed, rfSeed) {
  specs <- list()
  for (name in anchors) {
    if (name == "random_furthest") {
      for (j in seq_along(pivots)) {
        specs[[length(specs) + 1L]] <- list(
          label = paste0("random_furthest", j),
          s1    = as.integer(rfSeed(pivots[[j]]))
        )
      }
    } else {
      sd <- anchorSeed(name)
      specs[[length(specs) + 1L]] <- list(label = name, s1 = sd)
    }
  }
  if (length(specs) == 0L) {
    stop("no seed strategies to run; `pivots` is empty and the ensemble names ",
         "only `random_furthest`")
  }
  specs
}

#' Peripheral seed index for Gonzalez selection (distance matrix)
#'
#' @param d Square numeric distance matrix.
#' @param method Anchor name; see [MaxMinSeed()]. Also accepts `"first"` (1).
#' @return Integer seed index.
#' @keywords internal
.MaxMinSeed <- function(d, method) {
  switch(method,
    first   = 1L,
    centroid = stop("`centroid` seed requires coordinates; supply `points=` ",
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
    stop("Unknown seed method: ", method)
  )
}

#' Peripheral seed index for Gonzalez selection (coordinates)
#'
#' Coordinate counterpart of [.MaxMinSeed()]; each anchor is computed from the
#' `…FromPoints_cpp` primitives, bit-identical to the matrix path on Euclidean
#' data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param method Anchor name; see [MaxMinSeed()]. Also accepts `"first"` (1).
#' @return Integer seed index.
#' @keywords internal
.MaxMinSeedPoints <- function(points, method) {
  switch(method,
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
    stop("Unknown seed method: ", method)
  )
}

#' Peripheral seed index for Gonzalez farthest-first selection
#'
#' Returns the index of a single deterministic peripheral seed, the starting
#' point used by [FarFirst()] under the corresponding `seed` strategy. Useful
#' when composing a custom selection pass.
#'
#' Anchors:
#' \describe{
#'   \item{`"centroid"`}{The point farthest from the coordinate mean
#'     (`argmax ||x - x_bar||`), an `O(N * dim)` approximate diameter endpoint.
#'     Computed from coordinates, so it requires `points`; it is unavailable on
#'     the distance-matrix path, where `"peripheral"` serves the same role.}
#'   \item{`"peripheral"`}{Two sweeps: the point furthest from point 1, then
#'     the point furthest from that (a diameter-endpoint approximation), in
#'     `O(N)`. The only anchor reachable from a distance-column oracle (the
#'     function path of [FarFirst()]).}
#'   \item{`"random_furthest"`}{The point furthest from a random pivot, in
#'     `O(N)`. The pivot is drawn with the session RNG; set a seed
#'     (`set.seed()`) for a reproducible index.}
#'   \item{`"diameter"`}{A row endpoint of the diameter pair (the maximum
#'     pairwise distance). Degenerate data (`d_max <= 0`) falls back to 1.}
#'   \item{`"anti_medoid"`}{The point furthest from the 1-median (medoid).}
#'   \item{`"medoid"`}{The 1-median itself, `which.min(rowSums(d))`.}
#'   \item{`"rowsum"`}{The point maximising the sum of distances to all others
#'     (the 1-anti-median).}
#'   \item{`"rownorm"`}{The point maximising `sqrt(sum d^2)`, the L2 counterpart
#'     of `"rowsum"`.}
#' }
#'
#' @param d A `dist` object or square symmetric numeric matrix. Ignored when
#'   `points` is supplied.
#' @param points Optional `N x dim` numeric coordinate matrix; when supplied the
#'   seed is computed from coordinates in `O(N)` memory. Required for the
#'   `"centroid"` anchor, which has no distance-matrix form.
#' @param method Anchor name (see above). Default `"peripheral"`.
#' @return Integer seed index in `[1, N]`.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' MaxMinSeed(d, method = "diameter")
#' FarFirst(d, 5L, seed = MaxMinSeed(d, method = "diameter"))
#' @seealso [FarFirst()], which seeds and runs the greedy pass in one call.
#' @export
MaxMinSeed <- function(d = NULL, points = NULL,
                       method = c("peripheral", "centroid", "random_furthest",
                                  "diameter", "anti_medoid", "medoid", "rowsum",
                                  "rownorm")) {
  method <- match.arg(method)
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    return(.MaxMinSeedPoints(points, method))
  }
  d <- .AsDistMatrix(d)
  .MaxMinSeed(d, method)
}

#' Ensemble Gonzalez over cheap peripheral-anchor strategies (distance matrix)
#'
#' Runs Gonzalez from each requested peripheral anchor and returns the subset
#' maximising \eqn{T_k}. Internal driver for the ensemble path of [FarFirst()]
#' (triggered when `seed` is a character vector of length > 1). The
#' `"random_furthest"` token expands to one start per element of `pivots`. The
#' returned vector carries `strategy_results` and `winning_strategy` (character
#' vector of all tied-best strategies, with random starts labelled
#' `random_furthest1`, `random_furthest2`, ...) attributes.
#' @param d Square numeric distance matrix (already coerced).
#' @param n Integer subset size (`1 <= n < nrow(d)`).
#' @param anchors Character vector of anchor names.
#' @param pivots Integer vector of pivot indices the `"random_furthest"` token
#'   expands over (empty contributes none).
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsemble <- function(d, n, anchors = "peripheral", pivots = integer(0)) {
  d <- .AsDistMatrix(d)
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
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
  if (n >= nPts) return(seq_len(nPts))
  if (n < 1L)   return(integer(0))

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
      idx <- .MaximinFrom(d, n, first = s1)
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

  expanded <- .ExpandAnchors(anchors, pivots, AnchorSeed,
                             function(r) which.max(d[, r]))
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  strategyResults <- vector("list", length(expanded))
  names(strategyResults) <- labels
  for (i in seq_along(expanded)) {
    g <- RunGonz(expanded[[i]]$s1)
    strategyResults[[i]] <- list(
      s1  = expanded[[i]]$s1,
      idx = g$idx,
      t_k = g$tK
    )
  }

  bestI  <- 1L
  bestTk <- strategyResults[[1L]]$t_k
  for (i in seq_along(strategyResults)[-1L]) {
    tk <- strategyResults[[i]]$t_k
    if (is.na(tk)) next
    if (is.na(bestTk) || tk > bestTk) {
      bestI  <- i
      bestTk <- tk
    }
  }

  winners <- if (is.na(bestTk)) {
    bestI
  } else {
    which(vapply(strategyResults, function(r) isTRUE(r$t_k == bestTk),
                 logical(1L)))
  }

  result <- strategyResults[[bestI]]$idx
  attr(result, "strategy_results") <- strategyResults
  attr(result, "winning_strategy") <- labels[winners]
  result
}

#' Coordinate (matrix-free) multi-anchor Gonzalez ensemble
#'
#' Coordinate counterpart of [.GonzEnsemble()]; each anchor seed and the greedy
#' expansion are computed from `points` via the coordinate primitives, so the
#' returned indices and attributes match the matrix path on Euclidean data.
#' `pivots` indexes points directly, so the `"random_furthest"` starts also
#' match the matrix path.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param n Integer subset size.
#' @param anchors Character vector of anchor names.
#' @param pivots Integer vector of pivot indices the `"random_furthest"` token
#'   expands over (empty contributes none).
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsembleFromPoints <- function(points, n, anchors = .kDefaultEnsemble,
                                    pivots = integer(0)) {
  points <- .AsPointsMatrix(points)
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
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
  if (n >= nPts) return(seq_len(nPts))
  if (n == 0L)   return(integer(0))

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
      idx <- .MaximinFromPoints(points, n, first = s1)
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

  expanded <- .ExpandAnchors(
    anchors, pivots, AnchorSeed,
    function(r) which.max(EuclidColFromPoints_cpp(points, r))
  )
  labels   <- vapply(expanded, `[[`, character(1L), "label")
  strategyResults <- vector("list", length(expanded))
  names(strategyResults) <- labels
  for (i in seq_along(expanded)) {
    g <- RunGonz(expanded[[i]]$s1)
    strategyResults[[i]] <- list(
      s1  = expanded[[i]]$s1,
      idx = g$idx,
      t_k = g$tK
    )
  }

  bestI  <- 1L
  bestTk <- strategyResults[[1L]]$t_k
  for (i in seq_along(strategyResults)[-1L]) {
    tk <- strategyResults[[i]]$t_k
    if (is.na(tk)) next
    if (is.na(bestTk) || tk > bestTk) {
      bestI  <- i
      bestTk <- tk
    }
  }

  winners <- if (is.na(bestTk)) {
    bestI
  } else {
    which(vapply(strategyResults, function(r) isTRUE(r$t_k == bestTk),
                 logical(1L)))
  }

  result <- strategyResults[[bestI]]$idx
  attr(result, "strategy_results") <- strategyResults
  attr(result, "winning_strategy") <- labels[winners]
  result
}
