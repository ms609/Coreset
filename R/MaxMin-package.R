#' MaxMin: Maximum-Minimum Diversity and Dispersion Subset Selection
#'
#' Selects a representative subset of a fixed candidate set under one of three
#' distance-based objectives: the Max-Min Diversity Problem (MMDP, the discrete
#' *p*-dispersion objective: maximise the minimum pairwise distance within the
#' chosen subset); the Max-Mean Dispersion Problem (maximise the mean pairwise
#' distance, the subset size being free); and the discrete *k*-centre problem
#' (minimise the largest distance from any element to its nearest centre). The
#' solvers operate on a distance matrix, on Euclidean coordinates (without
#' materialising the matrix), or on an on-demand distance-column oracle.
#'
#' @section Max-Min diversity solvers:
#' \describe{
#'   \item{[FarFirst()]}{Greedy farthest-first selection from a distance matrix,
#'     a coordinate matrix, or a distance-column oracle (for spaces with no
#'     coordinate embedding), with a choice of peripheral seeding strategies and
#'     a robust ensemble default.}
#'   \item{[DropAdd()]}{DropAdd tabu search heuristic.}
#'   \item{[Grasp()]}{GRASP with path relinking.}
#'   \item{[ExactMaxMin()]}{Exact node-packing optimum (needs \pkg{highs}).}
#' }
#'
#' @section Max-Mean dispersion solver:
#' \describe{
#'   \item{[MaxMean()]}{Reinforcement-learning tabu search (RLTS).}
#' }
#'
#' @section k-centre solvers:
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
#' @section Relation to \pkg{maximin}:
#' Not to be confused with the CRAN package \pkg{maximin} (Sun & Gramacy),
#' which constructs continuous *space-filling designs* — it generates new
#' points in a coordinate region to maximise the minimum inter-point distance.
#' `MaxMin` instead *selects a subset* from a *fixed* candidate set under an
#' arbitrary distance, a combinatorial problem on a different footing.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @useDynLib MaxMin, .registration = TRUE
"_PACKAGE"
