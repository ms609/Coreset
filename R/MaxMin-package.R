#' MaxMin: Maximum-Minimum Diversity and Dispersion Subset Selection
#'
#' Selects a maximally dispersed subset of a fixed candidate set under the
#' Max-Min Diversity Problem (MMDP, the discrete *p*-dispersion objective):
#' maximise the minimum pairwise distance within the chosen subset. The solvers
#' operate on a distance matrix, on Euclidean coordinates (without
#' materialising the matrix), or on an on-demand distance-column oracle.
#'
#' @section Solvers:
#' \describe{
#'   \item{[Gonzalez()]}{Greedy farthest-first selection from a distance matrix,
#'     a coordinate matrix, or a distance-column oracle (for spaces with no
#'     coordinate embedding), with a choice of peripheral seeding strategies and
#'     a robust ensemble default.}
#'   \item{[DropAddTS()] / [DropAddTSPoints()]}{DropAdd tabu search heuristic.}
#'   \item{[ExactMaxMin()]}{Exact node-packing optimum (needs \pkg{highs}).}
#' #'   \item{[MinDist()]}{The k-centre objective (minimum pairwise distance).}
#' }
#'
#' @section Relation to \pkg{maximin}:
#' Not to be confused with the CRAN package \pkg{maximin} (Sun & Gramacy),
#' which constructs continuous *space-filling designs* — it generates new
#' points in a coordinate region to maximise the minimum inter-point distance.
#' `MaxMin` instead *selects a subset* from a *fixed* candidate set under an
#' arbitrary distance, a combinatorial problem on a different footing.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @useDynLib MaxMin, .registration = TRUE
"_PACKAGE"
