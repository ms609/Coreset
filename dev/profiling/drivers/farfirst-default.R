# Driver: FarFirst default coordinate path (seed = "random_furthest", 3 starts)
# Question: where does time go — seed-find / greedy pass / MinDist scoring — and
# how does the split move with n/N? Decides whether MaximinFrom_cpp is the right
# C++ entry point or whether the ensemble wrangling should move into C++.
#
# bare: see comment at foot (re-stamped each run)
.libdir <- tryCatch(readLines("dev/profiling/.libdir", n = 1L), error = function(e) NULL)
suppressMessages({
  if (!is.null(.libdir) && nzchar(.libdir)) library(MaxMin, lib.loc = .libdir) else library(MaxMin)
  library(bench)
})
set.seed(5813)

EuclidCol  <- MaxMin:::EuclidColFromPoints_cpp
Maximin    <- MaxMin:::MaximinFromPoints_cpp
MinPair    <- MaxMin:::.MinPairwiseFromPoints

`%then%` <- function(a, b) b  # readability noop

bench1 <- function(expr, iters = 25L) {
  e <- substitute(expr)
  env <- parent.frame()
  m <- bench::mark(eval(e, env), iterations = iters, check = FALSE,
                   filter_gc = FALSE)
  as.numeric(m$median) * 1000  # ms
}

run_case <- function(N, dim, n, npiv = 3L) {
  pts    <- matrix(rnorm(N * dim), ncol = dim)
  pivots <- sample.int(N, npiv)

  # Phase A — seed-find: which.max over one full Euclid column per pivot
  seedA <- function() {
    vapply(pivots, function(r) which.max(EuclidCol(pts, r)), integer(1L))
  }
  seeds <- seedA()

  # Phase B — greedy pass: one Gonzalez pass per seed
  passB <- function() {
    for (s in seeds) Maximin(pts, n, as.integer(s), 0L)
    invisible()
  }
  idxs <- lapply(seeds, function(s) Maximin(pts, n, as.integer(s), 0L))

  # Phase C — scoring: MinDist on each selection (k x k stats::dist)
  scoreC <- function() {
    for (idx in idxs) MinPair(pts, idx)
    invisible()
  }

  # Integrated default call (fixed pivots for determinism)
  full <- function() {
    FarFirst(n = n, points = pts, seed = "random_furthest", pivots = pivots)
  }

  a <- bench1(seedA())
  b <- bench1(passB())
  c <- bench1(scoreC())
  f <- bench1(full(), iters = 15L)

  data.frame(
    N = N, dim = dim, n = n, ratio = round(n / N, 3),
    seed_ms = round(a, 2), pass_ms = round(b, 2), score_ms = round(c, 2),
    full_ms = round(f, 2),
    score_pct_of_full = round(100 * c / f, 1),
    score_over_pass   = round(c / b, 3)
  )
}

N <- 6000L
for (dim in c(2L, 10L)) {
  cat(sprintf("\nFarFirst default breakdown  (N=%d, dim=%d, 3 starts)\n", N, dim))
  res <- do.call(rbind, lapply(
    c(60L, 300L, 1200L, 3000L),
    function(n) run_case(N, dim, n)
  ))
  print(res, row.names = FALSE)
}

cat("\nlegend: seed/pass/score are the 3 phases summed over the 3 starts;\n")
cat("full = end-to-end FarFirst; score_over_pass should track ~n/N.\n")
