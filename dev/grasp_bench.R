# Speedup of the compiled Grasp kernel vs the pure-R reference, at identical
# parameters and seed (so both do exactly the same work and return the same
# answer). Confirms bit-identity and reports the C++/R wall-clock ratio.
suppressPackageStartupMessages(library(Coreset))

bench <- function(n, m, mni, es, seed = 1L) {
  set.seed(123)
  pts <- matrix(rnorm(n * 3), ncol = 3)
  d <- as.matrix(dist(pts))

  set.seed(seed)
  tr <- system.time(ref <- Coreset:::.Grasp_R(d, m, max_no_improve = mni,
                                               elite_size = es))[["elapsed"]]
  set.seed(seed)
  tk <- system.time(ker <- Grasp(d, m, max_no_improve = mni,
                                   elite_size = es))[["elapsed"]]

  ok <- identical(ref$indices, ker$indices) &&
        identical(ref$objective, ker$objective) &&
        identical(ref$iters, ker$iters)
  cat(sprintf("n=%4d m=%3d mni=%3d es=%2d | iters=%4d | R=%6.2fs  C++=%6.3fs  speedup=%5.1fx | identical=%s\n",
              n, m, mni, es, ker$iters, tr, tk,
              if (tk > 0) tr / tk else NA_real_, ok))
  invisible(ok)
}

cat("=== Grasp: compiled kernel vs R reference (same seed/params) ===\n")
ok1 <- bench(100,  10, 50, 8)
ok2 <- bench(200,  20, 50, 8)
ok3 <- bench(300,  15, 40, 10)
stopifnot(ok1, ok2, ok3)
cat("\nAll runs bit-identical; speedups above.\n")
