# score.R

#' Minimum pairwise distance within a selection
#'
#' Returns the minimum pairwise distance among selected points \eqn{T_k}.
#' A set of points that is more dispersed will exhibit a higher value.
#'
#' @param d Pairwise distance matrix or `dist` object. Ignored when `points`
#'   is supplied.
#' @param idx Integer vector of selected row/col indices.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the score is computed from `stats::dist()` on the selected
#'   sub-coordinates only (`k x k`).
#' @return `MinDist()` returns a numeric specifying the minimum distance between
#' two selected points (or `NA_real_`, if `length(idx) < 2`).
#' @details The solvers in this package ([FarFirst()], [DropAdd()], [Grasp()])
#'  already attach the achieved \eqn{T_k} as a `score` attribute.
#'  `MinDist()` allows arbitrary selections to be scored, or an existing
#'  selection to be scored against a different distance matrix.
#' @seealso [FarFirst()], [DropAdd()], [Grasp()] and [ExactMaxMin()].
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' MinDist(d, FarFirst(5L, d))
#' @export
MinDist <- function(d = NULL, idx, points = NULL) {
  idx <- as.integer(idx)
  # NA in `idx` previously gave a silent NA on the matrix path but an error on
  # the coordinate path (F-605); duplicate indices made any selection score 0,
  # since the self-distance survives `diag(sub) <- Inf` (F-604). Reject both.
  if (anyNA(idx)) {
    stop("`idx` must not contain NA")
  }
  if (anyDuplicated(idx)) {
    stop("`idx` must not contain duplicate indices")
  }
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    return(.MinPairwiseFromPoints(points, idx))
  }
  d <- .AsDistMatrix(d)
  .SubsetScore(d, idx, objective = "min_pairwise")
}

#' Mean dispersion of a selection
#'
#' `MeanDist()` reports the sum of pairwise distances divided by the number of
#' selected elements,
#' \deqn{f(S) = \frac{\displaystyle\sum_{i < j,\, i,j \in S} d_{ij}}{|S|}},
#' the objective maximised by [MaxMean()].
#'
#' @param d Pairwise distance matrix or `dist` object.
#' @param idx Integer vector of selected row/col indices.
#' @return `MeanDist()` returns a numeric scalar, or `NA_real_` if
#'   `length(idx) < 2`.
#' @seealso [MaxMean()] which maximises this objective; [MinDist()] for the
#'   max-min (MMDP) analogue.
#' @examples
#' d <- dist(matrix(rnorm(20), ncol = 2))
#' selection <- MaxMean(d)
#' MeanDist(d, selection)
#' @export
MeanDist <- function(d, idx) {
  idx <- as.integer(idx)
  if (anyNA(idx)) stop("`idx` must not contain NA")
  if (anyDuplicated(idx)) stop("`idx` must not contain duplicate indices")
  d <- .AsDistMatrix(d)
  k <- length(idx)
  if (k < 2L) return(NA_real_)
  sub <- d[idx, idx, drop = FALSE]
  # Sum over unordered pairs, symmetrizing any asymmetry as the MaxMean()
  # C++ kernel does.
  # Keeps MeanDist() in agreement with MaxMean()'s reported score on the
  # asymmetric matrices .AsDistMatrix() silently accepts.
  # Return:
  (sum(sub) - sum(diag(sub))) / 2 / k
}
