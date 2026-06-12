# score.R

#' Minimum pairwise distance within a selection
#'
#' Returns the minimum pairwise distance among selected points \eqn{T_k}.
#' A set of points that is more dispersed will exhibit a higher value.
#'
#' @param idx Integer vector of selected row/col indices.
#' @param d Pairwise distance matrix or `dist` object. Ignored when `points`
#'   is supplied.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the score is computed from `stats::dist()` on the selected
#'   sub-coordinates only (`k x k`), never the full `N x N` matrix (`d` is then
#'   unused). For Euclidean data the result is identical to the matrix path.
#' @return `MinDist()` returns a numeric scalar; `NA_real_` if `length(idx) < 2`.
#' @details The solvers in this package ([FarFirst()], [DropAdd()], [Grasp()])
#'   already attach the achieved \eqn{T_k} as a `score` attribute, so
#'   `MinDist()` is mainly for scoring a selection produced elsewhere -- a
#'   matrix-free or externally generated index set -- or for re-scoring an
#'   existing selection against a different distance matrix.
#' @seealso [FarFirst()], [DropAdd()], [Grasp()] and [ExactMaxMin()], whose
#'   results already carry the objective.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' MinDist(FarFirst(5L, d), d)
#' @export
MinDist <- function(idx, d = NULL, points = NULL) {
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
