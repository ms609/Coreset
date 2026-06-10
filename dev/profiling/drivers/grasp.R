# Driver: Grasp refinement (dense-matrix only)
# Question: Where does Grasp time go — construction / local search (phase B)
# vs path relinking (phase C) — at small m and large m?
#
# Design:
#   n = 200  m in {10 (small), 50 (n/4), 100 (n/2)}  — large m stress case
#   n = 500  m = 10 only  (n/4 and n/2 take >>8 s: shrunk to stay <=8 s bare)
#   eliteSize varied at large m (3/5/8) to isolate PR contribution (PR ~ eliteSize^2)
#   plateau = 15L, seed = 1L throughout for determinism
#   Timing: median of 5 timed solver calls (matrix built once outside)
#   profvis: n=200, m=100, eliteSize=5, 3 reps (~7 s total)
#
# bare wall time (dry-run 2026-06-10, block-median):
#   n=200 m= 10 eliteSize=5 plateau=15 : ~  5.5 ms/call  (block of 20)
#   n=200 m= 50 eliteSize=5 plateau=15 : ~440   ms/call
#   n=200 m=100 eliteSize=3 plateau=15 : ~1010  ms/call
#   n=200 m=100 eliteSize=5 plateau=15 : ~2200  ms/call
#   n=200 m=100 eliteSize=8 plateau=15 : ~5450  ms/call  [within 8 s bound]
#   n=500 m= 10 eliteSize=5 plateau=15 : ~ 14.5 ms/call  (block of 20)
#   n=500 m=125 (n/4) eliteSize=5     : >>8 s  --> SHRUNK: omitted from large-m
#   n=500 m=250 (n/2) eliteSize=5     : >>8 s  --> SHRUNK: omitted from large-m

.libdir <- tryCatch(readLines("dev/profiling/.libdir", n = 1L), error = function(e) NULL)
suppressMessages({
  if (!is.null(.libdir) && nzchar(.libdir)) library(MaxMin, lib.loc = .libdir) else library(MaxMin)
})

# ---- helpers ---------------------------------------------------------------

med5 <- function(expr) {
  e <- substitute(expr)
  env <- parent.frame()
  times <- numeric(5L)
  for (i in seq_along(times)) {
    t0 <- proc.time()[[3L]]
    eval(e, env)
    times[i] <- proc.time()[[3L]] - t0
  }
  median(times) * 1000  # ms
}

run_case <- function(d, m, plateau, eliteSize, label, n_timing = 5L) {
  # warm-up
  set.seed(1)
  r <- Grasp(d, m = m, plateau = plateau, eliteSize = eliteSize)
  iters    <- attr(r, "iters")
  pr_calls <- attr(r, "pr_calls")
  score    <- attr(r, "score")

  times <- numeric(n_timing)
  for (i in seq_len(n_timing)) {
    t0 <- proc.time()[[3L]]
    set.seed(1)
    Grasp(d, m = m, plateau = plateau, eliteSize = eliteSize)
    times[i] <- proc.time()[[3L]] - t0
  }
  ms_med <- round(median(times) * 1000, 1)
  cat(sprintf("%-30s  n=%3d  m=%3d  eliteSize=%d  plateau=%2d  iters=%3d  pr_calls=%3d  ms=%7.1f  score=%.5f\n",
              label, nrow(d), m, eliteSize, plateau, iters, pr_calls, ms_med, score))
  invisible(list(n = nrow(d), m = m, eliteSize = eliteSize, plateau = plateau,
                 iters = iters, pr_calls = pr_calls, ms = ms_med, score = score))
}

# ---- build distance matrices once ------------------------------------------

set.seed(42)
n200 <- 200L
d200 <- matrix(runif(n200 * n200), n200, n200)
d200 <- (d200 + t(d200)) / 2; diag(d200) <- 0

n500 <- 500L
d500 <- matrix(runif(n500 * n500), n500, n500)
d500 <- (d500 + t(d500)) / 2; diag(d500) <- 0

# ---- timing table ----------------------------------------------------------

cat("\n=== Grasp timing table ===\n\n")
cat(sprintf("%-30s  %-5s  %-5s  %-10s  %-8s  %-7s  %-9s  %-10s  %-7s\n",
            "label", "n", "m", "eliteSize", "plateau", "iters", "pr_calls", "ms_median", "score"))
cat(strrep("-", 120), "\n")

results <- list()

# n=200, small m
results[[1]]  <- run_case(d200, m = 10L,  plateau = 15L, eliteSize = 5L, "n200_m10_es5")

# n=200, medium m (n/4)
results[[2]]  <- run_case(d200, m = 50L,  plateau = 15L, eliteSize = 5L, "n200_m50_es5")

# n=200, large m (n/2) — vary eliteSize to isolate PR
results[[3]]  <- run_case(d200, m = 100L, plateau = 15L, eliteSize = 3L, "n200_m100_es3")
results[[4]]  <- run_case(d200, m = 100L, plateau = 15L, eliteSize = 5L, "n200_m100_es5")
results[[5]]  <- run_case(d200, m = 100L, plateau = 15L, eliteSize = 8L, "n200_m100_es8")

# n=500, small m only  (large m >> 8 s, shrunk)
results[[6]]  <- run_case(d500, m = 10L,  plateau = 15L, eliteSize = 5L, "n500_m10_es5")

cat("\n")

# ---- scaling analysis ------------------------------------------------------

cat("=== Scaling analysis: m (small vs large) at fixed n=200, eliteSize=5 ===\n")
r_m10  <- results[[1]]; r_m50  <- results[[2]]; r_m100 <- results[[4]]

ratio_10_50   <- r_m50$ms  / r_m10$ms
ratio_10_100  <- r_m100$ms / r_m10$ms
ratio_50_100  <- r_m100$ms / r_m50$ms
# empirical exponent: log(ratio) / log(m_ratio)
exp_10_50  <- log(ratio_10_50)  / log(50 / 10)
exp_10_100 <- log(ratio_10_100) / log(100 / 10)
exp_50_100 <- log(ratio_50_100) / log(100 / 50)

cat(sprintf("  m=10 -> m=50  : %.1f x  (empirical m-exponent %.2f)\n", ratio_10_50,  exp_10_50))
cat(sprintf("  m=50 -> m=100 : %.1f x  (empirical m-exponent %.2f)\n", ratio_50_100, exp_50_100))
cat(sprintf("  m=10 -> m=100 : %.1f x  (empirical m-exponent %.2f)\n", ratio_10_100, exp_10_100))
cat("\n")

cat("=== PR isolation: eliteSize scaling at n=200, m=100, fixed plateau=15 ===\n")
r_es3 <- results[[3]]; r_es5 <- results[[4]]; r_es8 <- results[[5]]
# PR pairs: eliteSize*(eliteSize-1)
pr3 <- r_es3$pr_calls; pr5 <- r_es5$pr_calls; pr8 <- r_es8$pr_calls
dt_es3_es5 <- r_es5$ms - r_es3$ms  # extra time from (20-6)=14 more PR calls
dt_es5_es8 <- r_es8$ms - r_es5$ms  # extra time from (56-20)=36 more PR calls
cat(sprintf("  eliteSize 3->5: +%d PR calls => +%.0f ms  (%.1f ms/PR call)\n",
            pr5 - pr3, dt_es3_es5, dt_es3_es5 / (pr5 - pr3)))
cat(sprintf("  eliteSize 5->8: +%d PR calls => +%.0f ms  (%.1f ms/PR call)\n",
            pr8 - pr5, dt_es5_es8, dt_es5_es8 / (pr8 - pr5)))
cat(sprintf("  [check] ms_es5 = %.0f ms, with ~%.0f ms from PR and ~%.0f ms from constr+LS\n",
            r_es5$ms,
            r_es5$pr_calls * (dt_es3_es5 / (pr5 - pr3)),
            r_es5$ms - r_es5$pr_calls * (dt_es3_es5 / (pr5 - pr3))))
cat("\n")

# ---- phase discriminator: maxIter=0 isolates phase A + phase C ------------

cat("=== Phase discriminator: maxIter=0 at n=200, m=100, eliteSize=5 ===\n")
{
  phase_times <- numeric(5L)
  for (i in seq_len(5L)) {
    t0 <- proc.time()[[3L]]
    set.seed(1)
    Grasp(d200, m = 100L, plateau = 15L, eliteSize = 5L, maxIter = 0L)
    phase_times[i] <- proc.time()[[3L]] - t0
  }
  set.seed(1)
  r_noB <- Grasp(d200, m = 100L, plateau = 15L, eliteSize = 5L, maxIter = 0L)

  phaseA_only_times <- numeric(5L)
  for (i in seq_len(5L)) {
    t0 <- proc.time()[[3L]]
    set.seed(1)
    Grasp(d200, m = 100L, plateau = 15L, eliteSize = 1L, maxIter = 0L)
    phaseA_only_times[i] <- proc.time()[[3L]] - t0
  }

  ms_noB  <- round(median(phase_times) * 1000, 1)
  ms_phaseA1 <- round(median(phaseA_only_times) * 1000, 1)
  ms_full <- r_m100$ms

  phaseB_ms <- round(ms_full - ms_noB, 1)
  phaseC_ms <- round(ms_noB - 5 * ms_phaseA1, 1)  # ~5 builds in phase A
  phaseA_ms <- round(5 * ms_phaseA1, 1)

  cat(sprintf("  eliteSize=5 maxIter=0 (A+C)    : %6.1f ms  iters=%d pr_calls=%d\n",
              ms_noB, attr(r_noB, "iters"), attr(r_noB, "pr_calls")))
  cat(sprintf("  eliteSize=1 maxIter=0 (A only) : %6.1f ms/build\n", ms_phaseA1))
  cat(sprintf("  full run (A+B+C)               : %6.1f ms  iters=%d pr_calls=%d\n",
              ms_full, r_m100$iters, r_m100$pr_calls))
  cat(sprintf("\n  Phase A (~5 builds):   %5.1f ms  (%4.1f%% of full)\n",
              phaseA_ms, 100 * phaseA_ms / ms_full))
  cat(sprintf("  Phase B (GRASP iters): %5.1f ms  (%4.1f%%)\n",
              phaseB_ms, 100 * phaseB_ms / ms_full))
  cat(sprintf("  Phase C (PR calls):    %5.1f ms  (%4.1f%%)\n",
              phaseC_ms, 100 * phaseC_ms / ms_full))
}
cat("\n")

# ---- profvis triage --------------------------------------------------------

if (requireNamespace("profvis", quietly = TRUE)) {
  cat("=== profvis: n=200, m=100, eliteSize=5, plateau=15, 3 reps ===\n")
  p <- profvis::profvis({
    for (i in 1:3) {
      set.seed(1)
      Grasp(d200, m = 100L, plateau = 15L, eliteSize = 5L)
    }
  })
  out_html <- file.path("dev/profiling/drivers", "grasp-profvis.html")
  htmlwidgets::saveWidget(p, out_html, selfcontained = TRUE)
  cat("  Saved:", out_html, "\n\n")

  # Identify R vs C++ split by analysing call-stack depths.
  # profvis records each sample as a set of rows (one per stack frame) sharing
  # the same `time` value.  depth=1 is the outermost R frame (Grasp), depth=2
  # is the C++ kernel frame (Grasp_cpp).  A time-point with depth=2 labelled
  # "Grasp_cpp" means we are inside the .Call — counted as C++ time.
  # A time-point that only reaches depth=1 ("Grasp") means we are in the
  # R-wrapper overhead.
  pd <- p$x$message$prof
  if (!is.null(pd)) {
    all_times <- unique(pd$time)
    in_cpp <- vapply(all_times, function(t) {
      any(pd$label[pd$time == t] == "Grasp_cpp")
    }, logical(1L))
    cpp_pct <- 100 * mean(in_cpp)
    r_pct   <- 100 * mean(!in_cpp)
    cat(sprintf("  Total time-point samples  : %d\n", length(all_times)))
    cat(sprintf("  Inside C++ kernel         : %d (%.1f%%)\n", sum(in_cpp),  cpp_pct))
    cat(sprintf("  R-wrapper only            : %d (%.1f%%)\n", sum(!in_cpp), r_pct))
    # Report any R frame (depth=1, not Grasp) that appears in >2% of samples
    r_only_times <- all_times[!in_cpp]
    if (length(r_only_times) > 0) {
      r_labels <- pd$label[pd$time %in% r_only_times & pd$depth == 1L]
      freq <- sort(table(r_labels), decreasing = TRUE)
      cat("  R hotspots (>2% of total samples):\n")
      shown <- FALSE
      for (i in seq_len(length(freq))) {
        pct <- 100 * freq[i] / length(all_times)
        if (pct >= 2) {
          cat(sprintf("    %-40s %d (%.1f%%)\n", names(freq)[i], freq[i], pct))
          shown <- TRUE
        }
      }
      if (!shown) cat("    (none)\n")
    } else {
      cat("  R hotspots: none (0 samples outside C++ kernel)\n")
    }
  }
} else {
  cat("profvis not available — skipping HTML widget\n")
}

cat("\n=== Done ===\n")
