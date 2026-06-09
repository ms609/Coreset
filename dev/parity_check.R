# Bit-identical parity: MaxMin (new) vs FurthestPoint (old, installed).
suppressPackageStartupMessages({library(MaxMin); library(FurthestPoint)})
set.seed(42)
pts <- matrix(rnorm(80 * 4), ncol = 4)
d   <- dist(pts)
dm  <- as.matrix(d)
ok  <- TRUE
chk <- function(label, a, b) {
  pass <- identical(a, b)
  cat(sprintf("%-45s %s\n", label, if (pass) "OK" else "*** MISMATCH ***"))
  if (!pass) { ok <<- FALSE }
}
strip <- function(x) { attributes(x) <- NULL; x }

for (n in c(2L, 5L, 10L)) {
  # matrix path
  chk(sprintf("diameter   matrix n=%d", n),
      MaxMin::FarFirst(dm, n, seed = "diameter"),
      FurthestPoint::WideSampleGonzDiameter(dm, n))
  chk(sprintf("medoid     matrix n=%d", n),
      MaxMin::FarFirst(dm, n, seed = "medoid"),
      FurthestPoint::WideSampleMedoidFirst(dm, n))
  chk(sprintf("antimedoid matrix n=%d", n),
      MaxMin::FarFirst(dm, n, seed = "anti_medoid"),
      FurthestPoint::WideSampleAntiMedoid(dm, n))
  chk(sprintf("ensemble   matrix n=%d", n),
      strip(MaxMin::FarFirst(dm, n)),
      strip(FurthestPoint::WideSampleGonzEnsemble(dm, n)))
  chk(sprintf("ensemble-attr winner n=%d", n),
      attr(MaxMin::FarFirst(dm, n), "winning_strategy"),
      attr(FurthestPoint::WideSampleGonzEnsemble(dm, n), "winning_strategy"))
  chk(sprintf("first=1    matrix n=%d", n),
      MaxMin::FarFirst(dm, n, seed = 1L),
      FurthestPoint::FarFirst(dm, n, first = 1L))
  # points path
  chk(sprintf("diameter   points n=%d", n),
      MaxMin::FarFirst(n = n, points = pts, seed = "diameter"),
      FurthestPoint::WideSampleGonzDiameter(n = n, points = pts))
  chk(sprintf("medoid     points n=%d", n),
      MaxMin::FarFirst(n = n, points = pts, seed = "medoid"),
      FurthestPoint::WideSampleMedoidFirst(n = n, points = pts))
  chk(sprintf("antimedoid points n=%d", n),
      MaxMin::FarFirst(n = n, points = pts, seed = "anti_medoid"),
      FurthestPoint::WideSampleAntiMedoid(n = n, points = pts))
  chk(sprintf("ensemble   points n=%d", n),
      strip(MaxMin::FarFirst(n = n, points = pts)),
      strip(FurthestPoint::WideSampleGonzEnsemble(n = n, points = pts)))
}

# Distance-column oracle (now the function path of FarFirst())
colFn <- function(i) dm[, i]
chk("Gonzalez column-oracle first=1",
    MaxMin::FarFirst(colFn, 8L, N = nrow(dm), seed = 1L),
    FurthestPoint::GonzalezColumn(colFn, nrow(dm), 8L, first = 1L))
chk("Gonzalez column-oracle peripheral seed",
    MaxMin::FarFirst(colFn, 8L, N = nrow(dm)),
    FurthestPoint::GonzalezColumn(colFn, nrow(dm), 8L))

# DropAdd (deterministic, matrix + points)
chk("DropAdd matrix",
    MaxMin::DropAdd(dm, 8L, time_budget_s = 1, max_iter = 50L)$indices,
    FurthestPoint::DropAdd(dm, 8L, time_budget_s = 1, max_iter = 50L)$indices)
chk("DropAdd points",
    MaxMin::DropAdd(points = pts, m = 8L, timeBudgetS = 1, maxIter = 50L),
    FurthestPoint::DropAdd(points = pts, m = 8L, timeBudgetS = 1, maxIter = 50L))

# MinDist
chk("MinDist matrix", MaxMin::MinDist(dm, s0), FurthestPoint::MinDist(dm, s0))
chk("MinDist points", MaxMin::MinDist(idx = s0, points = pts),
    FurthestPoint::MinDist(idx = s0, points = pts))

# ExactMaxMin (needs highs)
if (requireNamespace("highs", quietly = TRUE)) {
  chk("ExactMaxMin",
      MaxMin::ExactMaxMin(dm[1:20, 1:20], 4L)$indices,
      FurthestPoint::ExactMaxMin(dm[1:20, 1:20], 4L)$indices)
} else cat("highs not installed; skipping ExactMaxMin\n")

# MaxMinSeed sanity vs internal anchors
chk("MaxMinSeed medoid matrix",
    MaxMin::MaxMinSeed(dm, method = "medoid"),
    as.integer(which.min(rowSums(dm))))

cat("\nOVERALL:", if (ok) "ALL PARITY CHECKS PASS" else "*** FAILURES ABOVE ***", "\n")
