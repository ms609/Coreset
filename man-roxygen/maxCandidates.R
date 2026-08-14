#' @param maxCandidates Integer: a composable-coreset tractability cap. When the
#'   number of candidate points `N` exceeds `maxCandidates`, [FarFirst()] thins
#'   the candidates to a `maxCandidates`-point coreset (with the deterministic,
#'   RNG-free `"peripheral"` seed, so no random stream is perturbed), the solver
#'   runs on the coreset, and the chosen indices are mapped back to the original
#'   numbering. This lets the solver produce a solution at scales where it would
#'   otherwise be intractable. `maxCandidates = 0` (or `Inf`) disables thinning
#'   and runs on the full problem; a cap at or above `N` is a no-op. A cap below
#'   `k` is an error. The default is <%= default %>, <%= default_basis %>.
#'   Thinning is **on by default**: an input larger than the cap is thinned (and
#'   a warning is emitted) unless `maxCandidates = 0` is passed.
