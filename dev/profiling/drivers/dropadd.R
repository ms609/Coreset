# Driver: DropAdd tabu-search profiling
# Round: DropAdd (area 2) — initial timing triage + construction/search decomposition
# Date:  2026-06-10
#
# bare wall time (single rep, 2026-06-10, Windows 10):
#   worst case n=20000 dim=10 m=10000: ~2.6 s/rep (construction ~1.9 s, search ~0.7 s)
#   matrix n=4000 m=2000: ~0.25 s/rep (construction ~0.16 s, search ~0.09 s)
#
# Key finding: construction dominates at large n/m.
#   Matrix O(n^2) seed (dropadd.cpp lines 50-58): 115 ms at n=4000 > 90 ms for 1500 iters.
#   Points O(m*n*dim) greedy: 1060-1880 ms at n=20000 m=10000 vs 390-680 ms search.
#   Search per-iter is FLAT in m (two O(n*dim) column passes; recompute branch negligible).
#
# Strategy:
#   - Fix iteration count (maxIter = 1500, plateau = .Machine$integer.max) so
#     every case does exactly 1500 tabu iterations.
#   - Also time with maxIter = 0L to isolate construction cost.
#   - Matrix path: n in {2000, 4000}, m in {10, round(n/2)}.
#   - Points path: n in {20000} dim in {2, 10}, m in {10, round(n/2)};
#                  also n=4000 dim=2 for direct cross-path comparison at same n.
#
# profvis triage (both include construction):
#   - Large-m matrix case  => dropadd-profvis-matrix-large.html
#   - Points n=20000 dim=2 case => dropadd-profvis-points-large.html

suppressMessages({
  library(Coreset)
  library(profvis)
  library(htmlwidgets)
})

set.seed(42)

MAX_ITER <- 1500L
PLATEAU  <- .Machine$integer.max

bench1 <- function(expr, iters = 8L) {
  e <- substitute(expr)
  env <- parent.frame()
  vals <- numeric(iters)
  for (i in seq_len(iters)) {
    t0 <- proc.time()[[3L]]
    eval(e, envir = env)
    vals[i] <- proc.time()[[3L]] - t0
  }
  median(vals) * 1000  # ms
}

cat("=== Building inputs (outside timed region) ===\n")

# Matrix inputs
pts2000 <- matrix(rnorm(2000 * 3), ncol = 3)
pts4000 <- matrix(rnorm(4000 * 3), ncol = 3)
cat("Building d2000 (n=2000)...")
d2000 <- as.matrix(dist(pts2000))
cat(" done\n")
cat("Building d4000 (n=4000)...")
d4000 <- as.matrix(dist(pts4000))
cat(" done\n")

# Points inputs
pts_4000_d2  <- matrix(rnorm(4000  * 2),  ncol = 2)
pts_20000_d2  <- matrix(rnorm(20000 * 2),  ncol = 2)
pts_20000_d10 <- matrix(rnorm(20000 * 10), ncol = 10)

cat("\n=== Timing cases (median of 8 reps, maxIter=1500 and maxIter=0) ===\n")

run_case <- function(path, dmat, pts, n_label, dim_label, m) {
  reps <- if (n_label >= 20000 && m >= 10000) 5L else 8L
  if (path == "matrix") {
    ms_tot <- bench1(DropAdd(d = dmat, m = m, maxIter = MAX_ITER,
                             plateau = PLATEAU), iters = reps)
    ms_con <- bench1(DropAdd(d = dmat, m = m, maxIter = 0L,
                             plateau = PLATEAU), iters = reps)
  } else {
    ms_tot <- bench1(DropAdd(points = pts, m = m, maxIter = MAX_ITER,
                             plateau = PLATEAU), iters = reps)
    ms_con <- bench1(DropAdd(points = pts, m = m, maxIter = 0L,
                             plateau = PLATEAU), iters = reps)
  }
  ms_srch <- ms_tot - ms_con
  data.frame(path = path, n = n_label, dim = dim_label, m = m,
             iters = MAX_ITER,
             ms_total   = round(ms_tot,  1),
             ms_construct = round(ms_con, 1),
             ms_search  = round(ms_srch, 1),
             ms_per_iter = round(ms_srch / MAX_ITER, 4),
             stringsAsFactors = FALSE)
}

results <- list()
cat("matrix n=2000, m=10 ...\n")
results[[1]]  <- run_case("matrix", d2000, NULL, 2000L, NA_integer_, 10L)
cat("matrix n=2000, m=1000 ...\n")
results[[2]]  <- run_case("matrix", d2000, NULL, 2000L, NA_integer_, 1000L)
cat("matrix n=4000, m=10 ...\n")
results[[3]]  <- run_case("matrix", d4000, NULL, 4000L, NA_integer_, 10L)
cat("matrix n=4000, m=2000 ...\n")
results[[4]]  <- run_case("matrix", d4000, NULL, 4000L, NA_integer_, 2000L)
cat("points n=4000, dim=2, m=10 ...\n")
results[[5]]  <- run_case("points", NULL, pts_4000_d2,   4000L,  2L, 10L)
cat("points n=4000, dim=2, m=2000 ...\n")
results[[6]]  <- run_case("points", NULL, pts_4000_d2,   4000L,  2L, 2000L)
cat("points n=20000, dim=2, m=10 ...\n")
results[[7]]  <- run_case("points", NULL, pts_20000_d2,  20000L, 2L, 10L)
cat("points n=20000, dim=2, m=10000 ...\n")
results[[8]]  <- run_case("points", NULL, pts_20000_d2,  20000L, 2L, 10000L)
cat("points n=20000, dim=10, m=10 ...\n")
results[[9]]  <- run_case("points", NULL, pts_20000_d10, 20000L, 10L, 10L)
cat("points n=20000, dim=10, m=10000 ...\n")
results[[10]] <- run_case("points", NULL, pts_20000_d10, 20000L, 10L, 10000L)

tab <- do.call(rbind, results)
cat("\n=== Timing table (ms_per_iter = search only, construction subtracted) ===\n")
print(tab, row.names = FALSE)

# Check for cases > 8s
over8 <- tab[tab$ms_total > 8000, ]
if (nrow(over8) > 0) {
  cat("\nWARNING: the following cases exceeded 8s:\n")
  print(over8[, c("path", "n", "dim", "m", "ms_total")], row.names = FALSE)
} else {
  cat("All cases within 8 s.\n")
}

# ===========================================================================
# profvis triage
# ===========================================================================
cat("\n=== profvis triage — matrix large-m case (n=4000, m=2000) ===\n")
K_prof <- 5L
p_matrix <- profvis::profvis({
  for (i in seq_len(K_prof)) {
    DropAdd(d = d4000, m = 2000L, maxIter = MAX_ITER, plateau = PLATEAU)
  }
})
out_mat <- file.path("dev/profiling/drivers", "dropadd-profvis-matrix-large.html")
htmlwidgets::saveWidget(p_matrix, out_mat, selfcontained = TRUE)
cat("Saved:", out_mat, "\n")

cat("\n=== profvis triage — points large-m case (n=20000, dim=2, m=10000) ===\n")
p_points <- profvis::profvis({
  for (i in seq_len(K_prof)) {
    DropAdd(points = pts_20000_d2, m = 10000L, maxIter = MAX_ITER,
            plateau = PLATEAU)
  }
})
out_pts <- file.path("dev/profiling/drivers", "dropadd-profvis-points-large.html")
htmlwidgets::saveWidget(p_points, out_pts, selfcontained = TRUE)
cat("Saved:", out_pts, "\n")

cat("\n=== Done ===\n")
