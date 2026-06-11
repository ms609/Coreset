# Does the NEW solver actually handle the four unseen larger targets (683-990)?
# This is the feasibility claim the whole optimisation rests on -- the dense
# solver could not build the packing matrix at these sizes. k=10 (Fig-1 regime)
# plus k=2,4,6 for the largest (vowel) to confirm the grid.
suppressMessages(library(MaxMin))
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")
big <- c("tc19_breastcancer", "tc13_pima", "tc17_vehicle", "tc18_vowel")
for (cn in big) {
  pts <- as.matrix(cases[[cn]][["points"]]); storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts)); N <- nrow(d)
  set.seed(1)
  t <- proc.time()[[3L]]
  r <- MaxMin::ExactMaxMin(d, m = 10L, timeBudgetS = 7200, progress = FALSE)
  cat(sprintf("%-18s n=%4d k=10  obj=%.6f  proven=%s  %.2fs\n",
              cn, N, r$objective, r$proven, proc.time()[[3L]] - t))
  flush.console()
}
