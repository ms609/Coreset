# compare-r2.R — baseline (old) vs new build for the round-2 changes.
old <- readRDS("dev/profiling/verify-r2-old.rds")
new <- readRDS("dev/profiling/verify-r2-new.rds")

keys <- union(names(old$corr), names(new$corr))
ident <- vapply(keys, function(k) identical(old$corr[[k]], new$corr[[k]]),
                logical(1L))
cat("=== Correctness: bit-identity vs baseline ===\n")
cat(sprintf("identical: %d / %d  (FarFirst + DropAdd + GraspPR)\n",
            sum(ident), length(keys)))
grp <- function(p) {
  ks <- keys[startsWith(keys, p)]
  sprintf("  %-9s %d/%d identical", p, sum(ident[ks]), length(ks))
}
cat(grp("ff"), "\n"); cat(grp("da_mx"), "\n"); cat(grp("da_pt"), "\n")
cat(grp("gp"), "\n")
bad <- keys[!ident]
if (length(bad)) { cat("MISMATCHES:\n"); for (b in head(bad, 20)) cat("   ", b, "\n") }

cat("\n=== Timing: old -> new ===\n")
for (k in names(new$timing)) {
  o <- old$timing[[k]]; nw <- new$timing[[k]]
  cat(sprintf("  %-34s %8.2f -> %8.2f ms  (%.2fx)\n",
              k, o$ms, nw$ms, o$ms / nw$ms))
}
