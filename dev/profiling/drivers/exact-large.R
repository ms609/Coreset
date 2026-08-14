# Stretch cells beyond the 56-cell ground-truth grid: the 683-990 targets at
# k = 10, plus the instances whose certifying probes the combinatorial
# reduction settles without an IP (spam n=4601 at every k; satellite n=6435
# at k=4). Each witness is re-validated against d.
# FP_CELLS overrides the cell list as comma-separated "case:k" pairs.
lib <- Sys.getenv("FP_LIB", "")
if (nzchar(lib)) library(MaxMin, lib.loc = lib) else library(MaxMin)
suppressMessages(library(highs))
options(MaxMin.progress = FALSE)
load("C:/Users/pjjg18/GitHub/furthest-point/data/cases.rda")

cells <- strsplit(strsplit(Sys.getenv(
  "FP_CELLS",
  paste("tc19_breastcancer:10,tc13_pima:10,tc17_vehicle:10,tc18_vowel:10,",
        "tc21_spam:4,tc21_spam:6,tc21_spam:10,tc14_satellite:4", sep = "")
), ",")[[1]], ":")

for (cell in cells) {
  cs <- cell[[1]]
  k <- as.integer(cell[[2]])
  pts <- as.matrix(cases[[cs]][["points"]])
  storage.mode(pts) <- "double"
  d <- as.matrix(stats::dist(pts))
  set.seed(4242)
  t0 <- proc.time()[[3L]]
  r <- ExactMaxMin(k, d, maxSeconds = as.numeric(Sys.getenv("FP_BUDGET", "600")))
  wall <- proc.time()[[3L]] - t0
  idx <- as.integer(r)
  sub <- d[idx, idx]
  diag(sub) <- Inf
  cat(sprintf("%-18s n=%4d k=%2d  score=%.6f proven=%s achieved=%.6f  %.1fs\n",
              cs, nrow(d), k, attr(r, "score"), attr(r, "proven"),
              min(sub), wall))
  rm(d, pts)
  invisible(gc(FALSE))
}
