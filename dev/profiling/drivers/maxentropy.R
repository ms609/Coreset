# Driver: MaxEntropy() maxdet selector -- step breakdown and end-to-end timing.
#
# The one O(n^3) cost is the kernel's PSD repair (R/maxentropy.R ->
# src/symeigen.cpp); everything else is O(n^2) or O(n^2 k). Runs against
# whichever Coreset is first on .libPaths(): set CORESET_LIB to A/B an
# installed build against the default library, e.g.
#   CORESET_LIB=dev/profiling/.vtune-lib-<stamp> Rscript dev/profiling/drivers/maxentropy.R
#
# Instance families (deterministic), all of which a user can hit:
#   euclid8  -- Euclidean, 8-D: numerically positive-definite kernel (the
#               Cholesky fast path; no eigendecomposition at all)
#   euclid2  -- Euclidean, 2-D: PSD in exact arithmetic but ~40% of the
#               computed eigenvalues are negative round-off (partial eigen,
#               negative side, p ~ 0.4 n)
#   pow12    -- Euclidean^1.2, 3-D: indefinite, ~97% negative (positive side,
#               p ~ 0.03 n)
#   cid      -- Clustering-information distance between random 40-leaf trees
#               (TreeDist): indefinite, ~14% negative at n = 1200; skipped if
#               TreeDist is not installed
# bare (new build, 2026-09-14): ~75 s; old build ~150 s.
lib <- Sys.getenv("CORESET_LIB", unset = NA)
if (!is.na(lib) && nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(Coreset))
cat("Coreset", as.character(packageVersion("Coreset")), "from", find.package("Coreset"), "\n")
set.seed(5813)
tm <- function(expr, R = 3L) {
  e <- parent.frame()
  median(replicate(R, system.time(eval(expr, e))[["elapsed"]]))
}
BenchDist <- function(n, dim, seed = 1L) {
  set.seed(seed)
  as.matrix(dist(matrix(rnorm(n * dim), ncol = dim)))
}
CidDist <- function(n, leaves = 40L) {
  if (!requireNamespace("TreeDist", quietly = TRUE)) return(NULL)
  set.seed(1)
  as.matrix(TreeDist::ClusteringInfoDistance(ape::rmtree(n, leaves)))
}
instances <- list(
  euclid8_500 = BenchDist(500L, 8L),
  euclid8_1200 = BenchDist(1200L, 8L),
  euclid2_1200 = BenchDist(1200L, 2L),
  pow12_1200 = BenchDist(1200L, 3L) ^ 1.2,
  cid_1200 = CidDist(1200L)
)
instances <- Filter(Negate(is.null), instances)

cat(sprintf("%-13s %5s %8s %8s %8s %8s %8s %8s\n", "instance", "n", "kernel", "prepare",
            "clip20", "cliphalf", "shift20", "trunc20"))
for (nm in names(instances)) {
  d <- instances[[nm]]
  n <- nrow(d)
  kern <- Coreset:::.MaxEntropyKernel(d)
  tKern <- tm(quote(Coreset:::.MaxEntropyKernel(d)))
  tPrep <- tm(quote(Coreset:::.MaxEntropyPrepare(kern, "clip")))
  kHalf <- n %/% 2L
  t20 <- tm(quote(MaxEntropy(20L, d)))
  tHalf <- tm(quote(MaxEntropy(kHalf, d)))
  tShift <- tm(quote(MaxEntropy(20L, d, repair = "shift")))
  tTrunc <- tm(quote(MaxEntropy(20L, d, repair = "truncate")))
  cat(sprintf("%-13s %5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", nm, n, tKern, tPrep,
              t20, tHalf, tShift, tTrunc))
}
