# coordinate-utils.R
#
# Matrix-free Euclidean primitives: a distance column and the diameter pair,
# computed from coordinates without materialising the N x N matrix. These back
# the coordinate paths of the solvers and are exposed so that callers can
# build further matrix-free Euclidean selectors on the same primitives, with
# results bit-identical to the distance-matrix path.

#' Euclidean distance column from coordinates (matrix-free)
#'
#' Returns the distances from point `i` to every point, computing them on the
#' fly from `points` (`O(N * dim)` time, `O(N)` memory) and reproducing
#' `stats::dist()`'s Euclidean bits exactly. This is column `i` of the
#' distance matrix, never materialising the matrix.
#'
#' @param points A numeric `N x dim` coordinate matrix (complete, no `NA`).
#' @param i Integer 1-based point index.
#' @return Numeric vector of length `N`; the self-distance at position `i` is 0.
#' @examples
#' pts <- matrix(c(0, 0, 3, 4), ncol = 2, byrow = TRUE)
#' PointColumn(pts, 1L)  # c(0, 5)
#' @seealso [PointDiameter()], [Gonzalez()].
#' @export
PointColumn <- function(points, i) {
  points <- .AsPointsMatrix(points)
  i <- as.integer(i)
  if (length(i) != 1L || is.na(i) || i < 1L || i > nrow(points)) {
    stop("`i` must be a single index in [1, nrow(points)]")
  }
  EuclidColFromPoints_cpp(points, i)
}

#' Diameter pair from coordinates (matrix-free)
#'
#' Finds the maximum pairwise Euclidean distance and one realising pair,
#' computed from `points` in `O(N^2 * dim)` time and `O(N)` memory. The row
#' index matches the first below-diagonal maximum under column-major
#' `which.max` on the distance matrix, so it agrees with the matrix path.
#'
#' @param points A numeric `N x dim` coordinate matrix (complete, no `NA`).
#' @return Numeric vector `c(d_max, i, j)`: the maximum distance and the 1-based
#'   indices of a pair achieving it.
#' @examples
#' pts <- matrix(c(0, 0, 3, 4, 1, 1), ncol = 2, byrow = TRUE)
#' PointDiameter(pts)
#' @seealso [PointColumn()], [Gonzalez()].
#' @export
PointDiameter <- function(points) {
  points <- .AsPointsMatrix(points)
  out <- DiameterFromPoints_cpp(points)
  stats::setNames(out, c("d_max", "i", "j"))
}
