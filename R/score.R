# score.R

#' Minimum pairwise distance within a selection (T_k = k-centre objective)
#'
#' Returns the minimum pairwise distance among selected points, the canonical
#' k-centre objective \eqn{T_k}. Higher values indicate a more spread-out
#' (better) selection.
#'
#' @param d Pairwise distance matrix or `dist` object. Ignored when `points`
#'   is supplied.
#' @param idx Integer vector of selected row/col indices.
#' @param points Optional `N x dim` numeric coordinate matrix. When supplied,
#'   the score is computed from `stats::dist()` on the selected
#'   sub-coordinates only (`k x k`), never the full `N x N` matrix (`d` is then
#'   unused). For Euclidean data the result is identical to the matrix path.
#' @return Numeric scalar; `NA_real_` if `length(idx) < 2`.
#' @examples
#' set.seed(1)
#' pts <- matrix(rnorm(60), ncol = 2)
#' d <- dist(pts)
#' TkScore(d, Gonzalez(d, 5L))
#' @export
TkScore <- function(d = NULL, idx, points = NULL) {
  if (!is.null(points)) {
    points <- .AsPointsMatrix(points)
    return(.MinPairwiseFromPoints(points, as.integer(idx)))
  }
  d <- .AsDistMatrix(d)
  .SubsetScore(d, idx, objective = "min_pairwise")
}
