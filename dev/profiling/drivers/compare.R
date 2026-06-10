# compare.R  — verdict from verify-old.rds vs verify-new.rds
old <- readRDS("dev/profiling/verify-old.rds")
new <- readRDS("dev/profiling/verify-new.rds")

keys <- intersect(names(old$corr), names(new$corr))
stopifnot(length(keys) == length(old$corr), length(keys) == length(new$corr))

# 1. NEW build: matrix-path vs points-path identical selection?
crossNew <- vapply(new$corr, function(c) isTRUE(c$cross_ident), logical(1L))
# 2. OLD build: same (baseline reference for the guarantee)
crossOld <- vapply(old$corr, function(c) isTRUE(c$cross_ident), logical(1L))
# 3. OLD vs NEW selection identity (did the optimisation change the answer?)
sameMat <- vapply(keys, function(k)
  identical(old$corr[[k]]$mat, new$corr[[k]]$mat), logical(1L))
samePts <- vapply(keys, function(k)
  identical(old$corr[[k]]$pts, new$corr[[k]]$pts), logical(1L))
# 4. free-T_k vs MinDist() reference error (NEW build)
tkErr <- vapply(new$corr, function(c)
  if (is.na(c$tkErr)) 0 else c$tkErr, numeric(1L))

cat("=== Correctness (", length(keys), "cases ) ===\n")
cat(sprintf("cross-path identical  : NEW %d/%d   OLD %d/%d\n",
            sum(crossNew), length(crossNew), sum(crossOld), length(crossOld)))
cat(sprintf("OLD==NEW selection    : matrix %d/%d   points %d/%d\n",
            sum(sameMat), length(keys), sum(samePts), length(keys)))
cat(sprintf("free-T_k vs MinDist    : max abs err = %.3e\n", max(tkErr)))

report_mismatch <- function(flag, lab) {
  bad <- names(flag)[!flag]
  if (length(bad)) {
    cat(sprintf("  %s mismatches (%d): \n", lab, length(bad)))
    for (b in head(bad, 12L)) cat("    ", b, "\n")
  }
}
report_mismatch(setNames(crossNew, names(new$corr)), "cross-path NEW")
report_mismatch(setNames(sameMat, keys), "OLD!=NEW matrix")
report_mismatch(setNames(samePts, keys), "OLD!=NEW points")

cat("\n=== Timing  (full FarFirst, 3 starts, N=6000) ===\n")
cat(sprintf("%4s %6s %5s | %18s | %18s\n", "dim", "n", "ratio",
            "points old->new ms", "matrix old->new ms"))
for (k in names(new$timing)) {
  o <- old$timing[[k]]; nw <- new$timing[[k]]
  spP <- o$points_ms / nw$points_ms
  spM <- o$matrix_ms / nw$matrix_ms
  cat(sprintf("%4d %6d %5.2f | %7.1f -> %6.1f (%4.2fx) | %7.1f -> %6.1f (%4.2fx)\n",
              nw$dim, nw$n, nw$ratio,
              o$points_ms, nw$points_ms, spP,
              o$matrix_ms, nw$matrix_ms, spM))
}
