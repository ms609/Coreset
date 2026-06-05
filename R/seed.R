# seed.R
#
# Peripheral seeding strategies for Gonzalez farthest-first selection, and the
# ensemble driver that runs several and keeps the best by TkScore(). The single
# anchors are exposed through MaxMinSeed(); Gonzalez(seed = ) selects among them
# (or the ensemble) and runs the greedy pass.

#' Peripheral seed index for Gonzalez selection (distance matrix)
#'
#' @param d Square numeric distance matrix.
#' @param method Anchor name; see [MaxMinSeed()]. Also accepts `"first"` (1).
#' @return Integer seed index.
#' @keywords internal
.MaxMinSeed <- function(d, method) {
  switch(method,
    first   = 1L,
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
      d_off <- d
      diag(d_off) <- -Inf
      d_max <- max(d_off)
      if (!is.finite(d_max) || d_max <= 0) {
        1L
      } else {
        as.integer(arrayInd(which.max(d_off), dim(d_off))[1L, 1L])
      }
    },
    peripheral = {
      s1 <- which.max(d[, 1L])
      as.integer(which.max(d[, s1]))
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
    stop("Unknown seed method: ", method)
  )
}

#' Peripheral seed index for Gonzalez farthest-first selection
#'
#' Returns the index of a single deterministic peripheral seed, the starting
#' point used by [Gonzalez()] under the corresponding `seed` strategy. Useful
#' when composing a custom selection pass.
#'
#' Anchors:
#' \describe{
#'   \item{`"peripheral"`}{Two sweeps: the point furthest from point 1, then
#'     the point furthest from that (a diameter-endpoint approximation). The
#'     only anchor reachable from a distance-column oracle ([GonzalezColumn()]).}
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
#'   seed is computed from coordinates in `O(N)` memory.
#' @param method Anchor name (see above). Default `"peripheral"`.
#' @return Integer seed index in `[1, N]`.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' MaxMinSeed(d, method = "diameter")
#' Gonzalez(d, 5L, first = MaxMinSeed(d, method = "diameter"))
#' @seealso [Gonzalez()], which seeds and runs the greedy pass in one call.
#' @export
MaxMinSeed <- function(d = NULL, points = NULL,
                       method = c("peripheral", "diameter", "anti_medoid",
                                  "medoid", "rowsum", "rownorm")) {
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
#' Runs Gonzalez from each requested deterministic peripheral anchor and returns
#' the subset maximising \eqn{T_k}. Internal driver for `Gonzalez(seed =
#' "ensemble")`; see [Gonzalez()] for the anchor descriptions. The returned
#' vector carries `strategy_results` and `winning_strategy` attributes.
#' @param d Square numeric distance matrix (already coerced).
#' @param n Integer subset size (`1 <= n < nrow(d)`).
#' @param anchors Character vector of anchor names.
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsemble <- function(d, n,
                          anchors = c("diameter", "anti_medoid",
                                      "rowsum", "rownorm")) {
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
    choices = c("diameter", "anti_medoid", "rowsum", "rownorm"),
    several.ok = TRUE
  ))
  nPts <- nrow(d)
  if (n >= nPts) return(seq_len(nPts))
  if (n == 0L)   return(integer(0))

  lazy <- new.env(parent = emptyenv())
  get_row_sums <- function() {
    lazy$row_sums <- lazy$row_sums %||% rowSums(d)
    lazy$row_sums
  }
  get_row_sq_sums <- function() {
    lazy$row_sq_sums <- lazy$row_sq_sums %||% rowSums(d ^ 2)
    lazy$row_sq_sums
  }
  get_d_offdiag <- function() {
    lazy$d_offdiag <- lazy$d_offdiag %||% { tmp <- d; diag(tmp) <- -Inf; tmp }
    lazy$d_offdiag
  }
  get_medoid <- function() {
    lazy$medoid <- lazy$medoid %||% as.integer(which.min(get_row_sums()))
    lazy$medoid
  }

  gonz_cache <- new.env(parent = emptyenv())
  run_gonz <- function(s1, mask_medoid) {
    key <- paste0(s1, if (mask_medoid) "M" else "U")
    gonz_cache[[key]] %||% {
      if (mask_medoid) {
        med   <- get_medoid()
        d_use <- d
        d_use[med, ] <- -Inf
        d_use[, med] <- -Inf
        diag(d_use) <- 0
        idx <- .MaximinFrom(d_use, n, first = s1)
      } else {
        idx <- .MaximinFrom(d, n, first = s1)
      }
      t_k <- if (length(idx) >= 2L) TkScore(d, idx) else NA_real_
      res  <- list(idx = idx, t_k = t_k)
      gonz_cache[[key]] <- res
      res
    }
  }

  anchor_seed <- function(name) {
    switch(name,
      diameter = {
        d_off <- get_d_offdiag()
        d_max <- max(d_off)
        if (!is.finite(d_max) || d_max <= 0) {
          list(s1 = 1L, mask = FALSE)
        } else {
          pair_ij <- arrayInd(which.max(d_off), dim(d_off))
          list(s1 = as.integer(pair_ij[1L, 1L]), mask = FALSE)
        }
      },
      anti_medoid = {
        med <- get_medoid()
        dist_from_medoid <- d[, med]
        dist_from_medoid[med] <- -Inf
        list(s1 = as.integer(which.max(dist_from_medoid)), mask = FALSE)
      },
      rowsum = list(s1 = as.integer(which.max(get_row_sums())), mask = FALSE),
      rownorm = list(s1 = as.integer(which.max(get_row_sq_sums())), mask = FALSE)
    )
  }

  strategy_results <- vector("list", length(anchors))
  names(strategy_results) <- anchors
  for (i in seq_along(anchors)) {
    seed <- anchor_seed(anchors[[i]])
    g    <- run_gonz(seed$s1, seed$mask)
    strategy_results[[i]] <- list(
      s1     = seed$s1,
      idx    = g$idx,
      t_k    = g$t_k,
      chosen = FALSE
    )
  }

  best_i  <- 1L
  best_tk <- strategy_results[[1L]]$t_k
  for (i in seq_along(strategy_results)[-1L]) {
    tk <- strategy_results[[i]]$t_k
    if (is.na(tk)) next
    if (is.na(best_tk) || tk > best_tk) {
      best_i  <- i
      best_tk <- tk
    }
  }
  strategy_results[[best_i]]$chosen <- TRUE

  result <- strategy_results[[best_i]]$idx
  attr(result, "strategy_results") <- strategy_results
  attr(result, "winning_strategy") <- anchors[[best_i]]
  result
}

#' Coordinate (matrix-free) four-anchor Gonzalez ensemble
#'
#' Coordinate counterpart of [.GonzEnsemble()]; each anchor seed and the greedy
#' expansion are computed from `points` via the coordinate primitives, so the
#' returned indices and attributes match the matrix path on Euclidean data.
#' @param points A `double` `N x dim` coordinate matrix.
#' @param n Integer subset size.
#' @param anchors Character vector of anchor names.
#' @return Integer vector of selected indices with attributes.
#' @keywords internal
.GonzEnsembleFromPoints <- function(points, n,
                                    anchors = c("diameter", "anti_medoid",
                                                "rowsum", "rownorm")) {
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
    choices = c("diameter", "anti_medoid", "rowsum", "rownorm"),
    several.ok = TRUE
  ))
  nPts <- nrow(points)
  if (n >= nPts) return(seq_len(nPts))
  if (n == 0L)   return(integer(0))

  lazy <- new.env(parent = emptyenv())
  get_row_sums <- function() {
    lazy$row_sums <- lazy$row_sums %||% RowSumsFromPoints_cpp(points)
    lazy$row_sums
  }
  get_row_sq_sums <- function() {
    lazy$row_sq_sums <- lazy$row_sq_sums %||% RowSqSumsFromPoints_cpp(points)
    lazy$row_sq_sums
  }
  get_diameter <- function() {
    lazy$diameter <- lazy$diameter %||% DiameterFromPoints_cpp(points)
    lazy$diameter
  }
  get_medoid <- function() {
    lazy$medoid <- lazy$medoid %||% as.integer(which.min(get_row_sums()))
    lazy$medoid
  }

  gonz_cache <- new.env(parent = emptyenv())
  run_gonz <- function(s1, mask_medoid) {
    key <- paste0(s1, if (mask_medoid) "M" else "U")
    gonz_cache[[key]] %||% {
      idx <- if (mask_medoid) {
        .MaximinFromPoints(points, n, first = s1, mask = get_medoid())
      } else {
        .MaximinFromPoints(points, n, first = s1)
      }
      t_k <- if (length(idx) >= 2L) TkScore(idx = idx, points = points) else NA_real_
      res <- list(idx = idx, t_k = t_k)
      gonz_cache[[key]] <- res
      res
    }
  }

  anchor_seed <- function(name) {
    switch(name,
      diameter = {
        diam <- get_diameter()
        if (!is.finite(diam[1L]) || diam[1L] <= 0) {
          list(s1 = 1L, mask = FALSE)
        } else {
          list(s1 = as.integer(diam[2L]), mask = FALSE)
        }
      },
      anti_medoid = {
        med <- get_medoid()
        dist_from_medoid <- EuclidColFromPoints_cpp(points, med)
        dist_from_medoid[med] <- -Inf
        list(s1 = as.integer(which.max(dist_from_medoid)), mask = FALSE)
      },
      rowsum = list(s1 = as.integer(which.max(get_row_sums())), mask = FALSE),
      rownorm = list(s1 = as.integer(which.max(get_row_sq_sums())), mask = FALSE)
    )
  }

  strategy_results <- vector("list", length(anchors))
  names(strategy_results) <- anchors
  for (i in seq_along(anchors)) {
    seed <- anchor_seed(anchors[[i]])
    g    <- run_gonz(seed$s1, seed$mask)
    strategy_results[[i]] <- list(
      s1     = seed$s1,
      idx    = g$idx,
      t_k    = g$t_k,
      chosen = FALSE
    )
  }

  best_i  <- 1L
  best_tk <- strategy_results[[1L]]$t_k
  for (i in seq_along(strategy_results)[-1L]) {
    tk <- strategy_results[[i]]$t_k
    if (is.na(tk)) next
    if (is.na(best_tk) || tk > best_tk) {
      best_i  <- i
      best_tk <- tk
    }
  }
  strategy_results[[best_i]]$chosen <- TRUE

  result <- strategy_results[[best_i]]$idx
  attr(result, "strategy_results") <- strategy_results
  attr(result, "winning_strategy") <- anchors[[best_i]]
  result
}
