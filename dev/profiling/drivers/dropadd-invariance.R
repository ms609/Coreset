# nCores-invariance check for the DropAdd points kernel
# (farfirst-invariance.R pattern): the same call must return the identical
# selection, score and secondary at every thread count. Sizes sit above
# DA_PAR_MIN (16384) so the parallel pass regions genuinely engage.
library(MaxMin)
cat("lib:", dirname(system.file(package = "MaxMin")), "\n")

Check <- function(label, th) {
  base <- NULL
  ok <- TRUE
  for (nc in c(1L, 2L, 8L)) {
    options(mc.cores = nc)
    r <- th()
    v <- list(idx = as.integer(r), score = as.numeric(attr(r, "score")),
              sec = as.numeric(attr(r, "secondary")),
              iters = as.integer(attr(r, "iters")))
    if (is.null(base)) base <- v
    else if (!identical(base, v)) ok <- FALSE
  }
  options(mc.cores = NULL)
  cat(sprintf("%-46s nCores {1,2,8}: %s\n", label,
              if (ok) "INVARIANT" else "DIVERGED"))
  ok
}

set.seed(42)
big10 <- matrix(rnorm(20000L * 10L), ncol = 10L)
big2  <- matrix(rnorm(20000L * 2L),  ncol = 2L)

ok <- all(
  Check("points n=2e4 d10 k=10 plateau=200", function()
    DropAdd(10L, points = big10, plateau = 200L, maxCandidates = 0L)),
  Check("points n=2e4 d2 k=50 plateau=200", function()
    DropAdd(50L, points = big2, plateau = 200L, maxCandidates = 0L)),
  Check("points n=2e4 d10 k=2000 plateau=60", function()
    DropAdd(2000L, points = big10, plateau = 60L, maxCandidates = 0L))
)
cat(if (ok) "RESULT: NCORES-INVARIANT\n" else "RESULT: DIVERGED\n")
if (!ok) quit(status = 1L)
