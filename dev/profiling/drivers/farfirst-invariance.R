# nCores-invariance check for the FarFirst kernels (grasp-invariance.R
# pattern): the same call must return identical selections and scores at
# every thread count. Sizes chosen to engage the greedy pass's parallel path
# (N >= 32768) and the O(N^2) anchor primitives.
library(MaxMin)
cat("lib:", dirname(system.file(package = "MaxMin")), "\n")

Check <- function(label, th) {
  base <- NULL
  ok <- TRUE
  for (nc in c(1L, 2L, 8L)) {
    options(mc.cores = nc)
    r <- th()
    v <- list(idx = as.integer(r), score = as.numeric(attr(r, "score")))
    if (is.null(base)) base <- v
    else if (!identical(base, v)) ok <- FALSE
  }
  options(mc.cores = NULL)
  cat(sprintf("%-46s nCores {1,2,8}: %s\n", label,
              if (ok) "INVARIANT" else "DIVERGED"))
  ok
}

set.seed(99)
big2  <- matrix(rnorm(50000L * 2L),  ncol = 2L)
big10 <- matrix(rnorm(50000L * 10L), ncol = 10L)
mid   <- matrix(rnorm(1200L * 5L),  ncol = 5L)
dmid  <- as.matrix(dist(mid))

ok <- all(
  Check("points N=5e4 dim=2 k=400 (parallel pass)", function()
    FarFirst(400L, points = big2, strategy = 5L)),
  Check("points N=5e4 dim=10 k=400 (parallel pass)", function()
    FarFirst(400L, points = big10, strategy = 5L)),
  Check("matrix N=1200 anchors ensemble", function()
    FarFirst(60L, dmid, strategy = c("diameter", "anti_medoid", "rownorm"))),
  Check("points N=1200 anchors ensemble", function()
    FarFirst(60L, points = mid,
             strategy = c("diameter", "anti_medoid", "medoid", "rowsum",
                          "rownorm", "peripheral"))),
  Check("points N=5e4 default ensemble (seeded)", function() {
    set.seed(7)
    FarFirst(120L, points = big2)
  })
)
cat(if (ok) "RESULT: NCORES-INVARIANT\n" else "RESULT: DIVERGED\n")
