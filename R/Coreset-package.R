#' Coreset: Discrete Diversity, Dispersion, and Coverage Subset Selection
#'
#' Selects a representative subset of a fixed candidate set.
#'
#' @section Max-Min diversity solvers:
#' The Max-Min Diversity Problem (MMDP) maximises the minimum pairwise distance
#'  within a subset (the discrete *p*-dispersion objective).
#' \describe{
#'   \item{[FarFirst()]}{Greedy farthest-first selection from a distance matrix,
#'     a coordinate matrix, or a distance-column oracle (for spaces with no
#'     coordinate embedding), with a choice of peripheral seeding strategies and
#'     a robust ensemble default.}
#'   \item{[DropAdd()]}{DropAdd tabu search heuristic.}
#'   \item{[Grasp()]}{GRASP with path relinking.}
#'   \item{[ExactMaxMin()]}{Exact node-packing optimum.}
#' }
#'
#' @section Max-Mean dispersion solver:
#' The Max-Mean Dispersion Problem selects a subset of a size that maximises
#' the mean pairwise distance.
#' \describe{
#'   \item{[MaxMean()]}{Reinforcement-learning tabu search.}
#' }
#'
#' @section k-centre solvers:
#' The discrete *k*-centre problem minimises the largest distance from any
#' element to its nearest selected element ('centre').
#' \describe{
#'   \item{[KCentre()]}{CDSh covering heuristic.}
#'   \item{[ExactKCentre()]}{Exact minimum-cover optimum (needs \pkg{highs}).}
#' }
#'
#' @section Scoring:
#' \describe{
#'   \item{[MinDist()]}{Minimum pairwise distance (the max-min objective).}
#'   \item{[MeanDist()]}{Mean pairwise dispersion (the max-mean objective).}
#'   \item{[KCentreRadius()]}{Covering radius (the k-centre objective).}
#' }
#'
#' @section Options:
#' The solvers need `d[i, j]` and `d[j, i]` to agree exactly, being free to
#' read whichever is the cheaper memory access. A matrix that misses this only
#' through rounding has its triangles averaged, with a warning;
#' `options(Coreset.symmetryTolerance = )` sets how large a discrepancy — scaled
#' by `max(1, |d[i, j]|, |d[j, i]|)` — is repaired rather than refused. It
#' defaults to `100 * .Machine$double.eps`, as R's own `isSymmetric()` does;
#' set it to `0` to have any inexact matrix refused.
#'
#' @section Relation to \pkg{maximin}:
#' Not to be confused with the CRAN package \pkg{maximin}, which constructs
#' continuous space-filling designs by generating *new* points in a coordinate
#' region to maximise the minimum inter-point distance.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @useDynLib Coreset, .registration = TRUE
"_PACKAGE"
