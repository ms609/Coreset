# Capture the OLD (installed) ExactMaxMin proven optima + timings as the
# regression baseline, before installing the new build over it.
suppressMessages(library(MaxMin))
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
cs14 <- c("tc20_zoo","tc7_ring","tc9_iris","tc8_highdim_gaussians","tc3_multiscale",
          "tc1_uniform","tc2_two_unequal","tc6_density_gradient","tc16_sonar",
          "tc5_clusters_outliers","tc10_glass","tc4_hierarchical",
          "tc22_penguins","tc11_ionosphere")
ks <- c(2L, 4L, 6L, 10L)
rows <- list()
for (cs in cs14) {
  pts <- as.matrix(cases[[cs]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts))
  for (k in ks) {
    t <- proc.time()[[3L]]
    r <- MaxMin::ExactMaxMin(d, k, maxSeconds = 600)
    rows[[length(rows)+1L]] <- data.frame(case = cs, n = nrow(d), k = k,
      objective = r$objective, proven = r$proven,
      time_s = proc.time()[[3L]] - t, stringsAsFactors = FALSE)
    cat(sprintf("OLD %-16s k=%2d  obj=%.6f  %.2fs\n", cs, k, r$objective,
                rows[[length(rows)]]$time_s))
  }
}
old <- do.call(rbind, rows)
saveRDS(old, "dev/profiling/oldopt.rds")
cat(sprintf("\nsaved %d rows; total OLD solve time %.1fs\n", nrow(old), sum(old$time_s)))
