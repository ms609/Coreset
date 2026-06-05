# Shared setup for the MaxMin performance-regression benchmarks.
#
# Each bench-*.R sources this file, builds a realistic fixed instance (outside
# the timed expression), and calls Benchmark(<call>). The harness saves a
# <call>.bench.Rds for benchmark/_compare_results.R to diff PR vs target.
#
# All refinement methods are driven by their DETERMINISTIC stopping rule
# (max_no_improve), never a wall-clock budget: with a time budget every run
# would measure the budget, not the work, so a regression would be invisible.

library("MaxMin")

`%||%` <- function(x, y) if (is.null(x)) y else x

Benchmark <- function(..., min_iterations = NULL, min_time = NULL) {
  # Pass ... straight to bench::mark to preserve the captured expression
  # (do.call would evaluate it first and lose the label).
  result <- if (is.null(min_time)) {
    bench::mark(..., min_iterations = min_iterations %||% 3, time_unit = "us",
                check = FALSE)
  } else {
    bench::mark(..., min_iterations = min_iterations %||% 3,
                min_time = min_time, time_unit = "us", check = FALSE)
  }
  if (interactive()) {
    print(result)
  } else {
    fileroot <- gsub("[\"']", "",
                     gsub("[\\(\\):, /]", "_", as.character(result$expression)))
    .FileName <- function(fileRoot, i) {
      paste0(c(fileroot, i, "bench.Rds"), collapse = ".")
    }
    i <- double(0)
    while (file.exists(.FileName(fileroot, i))) {
      if (length(i) == 0) i <- 0
      i <- 1 + i
    }
    saveRDS(result, .FileName(fileroot, i))
  }
}

# Deterministic Euclidean point cloud, n x dim.
BenchPoints <- function(n, dim, seed = 1L) {
  set.seed(seed)
  matrix(stats::rnorm(n * dim), ncol = dim)
}

# Deterministic dense distance matrix from a point cloud.
BenchDist <- function(n, dim, seed = 1L) {
  as.matrix(stats::dist(BenchPoints(n, dim, seed)))
}
