# dev/bench-dropadd-oracle.R
#
# Scratch benchmark (NOT part of the test suite) for DropAdd()'s distance-column
# oracle path, using the motivating consumer's metric: tree-to-tree distances
# from TreeDist.
#
# Question: does a *column* oracle recover the efficiency of a bulk all-pairs
# routine, or is its value only that it never materialises the N x N matrix?
#
# Compares, at a couple of representative (N, k):
#   (a) bulk all-pairs call + DropAdd matrix path (the status quo);
#   (b) oracle path, "good closure"  -- splits precomputed once, looked up;
#   (c) oracle path, "good closure" + memoized columns;
#   (d) oracle path, "bad closure"   -- splits re-derived on every call
#       (measured per column only; a full run at this cost is not worth the
#       wall-clock).
#
# Run against an INSTALLED build in a fresh session -- pkgload::load_all()
# compiles the C++ matrix kernels at -O0 and would flatter the pure-R oracle
# path (AGENTS.md).
#
#   R --vanilla -f dev/bench-dropadd-oracle.R
#
# Results and interpretation: dev/dropadd-oracle-treedist.md

suppressPackageStartupMessages({
  library("MaxMin")
  library("TreeDist")
  library("TreeTools")
})

stopifnot(packageVersion("MaxMin") >= "0.0.0.9005")
options(MaxMin.progress = FALSE)

kPlateau <- 200L      # keep the runs finishable; see the note on plateau below
kNTips   <- 60L

# A counting wrapper, so "oracle calls" is measured, not asserted.
CountingOracle <- function(fn) {
  n <- new.env(parent = emptyenv())
  n$calls <- 0L
  list(
    fn = function(i) {
      n$calls <- n$calls + 1L
      fn(i)
    },
    Calls = function() n$calls
  )
}

# --- the three closure shapes ----------------------------------------------

# GOOD: derive each tree's splits ONCE, in the enclosing environment; each call
# then only pays for N already-extracted comparisons.
GoodClosure <- function(trees) {
  splits <- as.Splits(trees)
  function(i) as.numeric(ClusteringInfoDistance(splits[[i]], splits))
}

# GOOD + MEMO: as above, but each column is computed at most once, capping the
# oracle at N distinct evaluations at the cost of an N x (touched) cache.
MemoClosure <- function(trees) {
  splits <- as.Splits(trees)
  cache <- new.env(parent = emptyenv())
  function(i) {
    key <- as.character(i)
    if (is.null(cache[[key]])) {
      cache[[key]] <- as.numeric(ClusteringInfoDistance(splits[[i]], splits))
    }
    cache[[key]]
  }
}

# BAD: hand the raw trees over every time, so both sides are re-parsed per call.
BadClosure <- function(trees) {
  function(i) as.numeric(ClusteringInfoDistance(trees[[i]], trees))
}

Timed <- function(expr) {
  t <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, seconds = t)
}

Row <- function(...) {
  cat(sprintf("%-26s %10s %12s %14s %12s\n", ...))
}

for (N in c(200L, 400L)) {
  set.seed(1)
  trees <- lapply(seq_len(N), function(i) RandomTree(kNTips, root = FALSE))
  class(trees) <- "multiPhylo"

  for (k in c(10L, 25L)) {
    cat(sprintf("\n===== N = %d trees (%d tips), k = %d, plateau = %d =====\n",
                N, kNTips, k, kPlateau))
    Row("path", "seconds", "oracle_calls", "pair_evals", "peak_MB")

    # (a) bulk all-pairs + matrix path -------------------------------------
    bulk <- Timed(as.matrix(ClusteringInfoDistance(trees)))
    dmat <- bulk$value
    solve <- Timed(DropAdd(k, dmat, plateau = kPlateau, maxCandidates = 0L))
    Row("bulk matrix", sprintf("%.1f", bulk$seconds + solve$seconds), "-",
        format(N * (N - 1L) / 2L, big.mark = ","),
        sprintf("%.1f", N * N * 8 / 2^20))
    cat(sprintf("%-26s %10s (build %.1f s + search %.2f s), score %.4f\n", "",
                "", bulk$seconds, solve$seconds, attr(solve$value, "score")))

    # (b) oracle, good closure ---------------------------------------------
    oc <- CountingOracle(GoodClosure(trees))
    good <- Timed(DropAdd(k, oc$fn, N = N, plateau = kPlateau))
    Row("oracle (good closure)", sprintf("%.1f", good$seconds),
        format(oc$Calls(), big.mark = ","),
        format(oc$Calls() * N, big.mark = ","),
        sprintf("%.1f", N * 8 * 5 / 2^20))
    cat(sprintf("%-26s %10s score %.4f, same subset as matrix path: %s\n", "", "",
                attr(good$value, "score"),
                identical(as.integer(good$value), as.integer(solve$value))))

    # (c) oracle, good closure + memoisation --------------------------------
    memoFn <- MemoClosure(trees)
    distinct <- new.env(parent = emptyenv())
    distinct$seen <- integer(0)
    counted <- function(i) {
      distinct$seen <- c(distinct$seen, i)
      memoFn(i)
    }
    memo <- Timed(DropAdd(k, counted, N = N, plateau = kPlateau))
    nDistinct <- length(unique(distinct$seen))
    Row("oracle (memoized)", sprintf("%.1f", memo$seconds),
        sprintf("%s (%s distinct)", format(length(distinct$seen), big.mark = ","),
                nDistinct),
        format(nDistinct * N, big.mark = ","),
        sprintf("%.1f", nDistinct * N * 8 / 2^20))

    # (d) bad closure, per-column cost only ---------------------------------
    badFn <- BadClosure(trees)
    perBad <- Timed(for (i in 1:5) badFn(i))$seconds / 5
    goodFn <- GoodClosure(trees)
    perGood <- Timed(for (i in 1:5) goodFn(i))$seconds / 5
    cat(sprintf("%-26s per column: good %.3f s, bad %.3f s, bulk-amortised %.3f s\n",
                "per-call costs", perGood, perBad, bulk$seconds / N))
  }
}

cat("\nDone.\n")
