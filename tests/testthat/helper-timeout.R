# Defensive backstop for tests that exercise a wall-clock stopping rule.
#
# These tests deliberately disable every *deterministic* stopping criterion
# (plateau, maxIter set to .Machine$integer.max) so that only the time budget
# can end the search. If the budget logic ever regresses into a non-terminating
# loop, the bare test would hang until the CI job's own (multi-hour) timeout.
# Wrapping the call in a tight elapsed limit turns that into a fast, legible
# failure instead.
#
# setTimeLimit() fires at R evaluation checkpoints; the compiled kernels also
# poll Rcpp::checkUserInterrupt(), which gives the limit a checkpoint inside the
# C++ loops too. Any R-level reference path (e.g. .Grasp_R) raises the limit as
# an ordinary error, which is what we catch here.
#
# NB: only safe because the budget routines no longer call setTimeLimit() with a
# non-positive `elapsed` -- that would *remove* this outer limit (R semantics)
# and defeat the backstop.
expect_returns_within <- function(expr, limit = 5) {
  setTimeLimit(elapsed = limit, transient = TRUE)
  on.exit(setTimeLimit())
  tryCatch(
    expr,
    error = function(e) {
      if (grepl("elapsed time limit|time limit", conditionMessage(e),
                ignore.case = TRUE)) {
        stop(sprintf(
          "call did not return within %g s -- likely a non-terminating loop",
          limit), call. = FALSE)
      }
      stop(e)
    }
  )
}
